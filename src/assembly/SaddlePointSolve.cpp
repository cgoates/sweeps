// maxwell_eig.cpp
// Minimal PETSc / SLEPc C++ skeleton implementing matrix-free projected operators
// Build: see notes above. Fill the "ASSEMBLE" sections.

#include <petscksp.h>
#include <slepceps.h>
#include <iostream>
#include <memory>
#include <vector>
#include <SLEPcUtils.hpp>

// Context passed to shell Mat mult callbacks
struct ShellCtx
{
    Mat A;          // assembled A (edge bilinear)
    Mat M;          // assembled edge mass
    Mat B;          // coupling B (edge->node)
    Mat S_explicit; // optional: S = B * M^{-1} * B^T (may be MPI-distributed)
    bool S_is_explicit = false;

    // KSP solvers to apply M^{-1} and S^{-1}
    KSP ksp_M; // solves M * y = rhs (edge)
    KSP ksp_S; // solves S * x = rhs (nodal) — may act on S_explicit or implicit application

    // Temporary Vecs reused inside apply_P (created on same communicator as edge vector)
    Vec tmp_Bv;   // B * v  (nodal)
    Vec tmp_Sx;   // x solved from S x = Bv (nodal)
    Vec tmp_Bt_x; // B^T * x (edge)
    Vec tmp_My2;  // result of M^{-1} * (B^T x) (edge)
};

// Helper: apply projector P = I - M^{-1} B^T S^{-1} B
// Uses ShellCtx stored as MatShell context user pointer
PetscErrorCode Apply_P( Mat shellMat, Vec v, Vec out )
{
    PetscFunctionBeginUser;
    void* ctx_void = nullptr;
    MatShellGetContext( shellMat, &ctx_void );
    ShellCtx* ctx = static_cast<ShellCtx*>( ctx_void );
    PetscErrorCode ierr;

    // out = v (copy)
    ierr = VecCopy( v, out );
    CHKERRQ( ierr );

    // tmp_Bv = B * v
    ierr = MatMult( ctx->B, v, ctx->tmp_Bv );
    CHKERRQ( ierr );

    // Solve S * tmp_Sx = tmp_Bv  (ksp_S)
    ierr = KSPSolve( ctx->ksp_S, ctx->tmp_Bv, ctx->tmp_Sx );
    CHKERRQ( ierr );

    // tmp_Bt_x = B^T * tmp_Sx
    ierr = MatMultTranspose( ctx->B, ctx->tmp_Sx, ctx->tmp_Bt_x );
    CHKERRQ( ierr );

    // tmp_My2 = solve M * tmp_My2 = tmp_Bt_x
    ierr = KSPSolve( ctx->ksp_M, ctx->tmp_Bt_x, ctx->tmp_My2 );
    CHKERRQ( ierr );

    // out = out - tmp_My2
    ierr = VecAXPY( out, -1.0, ctx->tmp_My2 );
    CHKERRQ( ierr );

    PetscFunctionReturn( 0 );
}

PetscErrorCode Apply_S(Mat shell_S, Vec x, Vec y)
{
    PetscFunctionBeginUser;
    void* ctx_void;
    MatShellGetContext(shell_S, &ctx_void);
    ShellCtx* ctx = static_cast<ShellCtx*>(ctx_void);
    PetscErrorCode ierr;
    
    // Need temporary vectors
    Vec temp_edge;  // Size 880 (edge space)
    MatCreateVecs(ctx->B, &temp_edge, NULL);
    
    // y = B * M^(-1) * B^T * x
    // Step 1: temp_edge = B^T * x
    ierr = MatMultTranspose(ctx->B, x, temp_edge); CHKERRQ(ierr);
    
    // Step 2: solve M * temp_edge2 = temp_edge  (i.e., temp_edge2 = M^(-1) * temp_edge)
    Vec temp_edge2;
    VecDuplicate(temp_edge, &temp_edge2);
    ierr = KSPSolve(ctx->ksp_M, temp_edge, temp_edge2); CHKERRQ(ierr);
    
    // Step 3: y = B * temp_edge2
    ierr = MatMult(ctx->B, temp_edge2, y); CHKERRQ(ierr);
    
    VecDestroy(&temp_edge);
    VecDestroy(&temp_edge2);
    PetscFunctionReturn(0);
}

// y <- P^T * A * P * x
PetscErrorCode AP_mult( Mat shellMat, Vec x, Vec y )
{
    PetscFunctionBeginUser;
    void* ctx_void = nullptr;
    MatShellGetContext( shellMat, &ctx_void );
    ShellCtx* ctx = static_cast<ShellCtx*>( ctx_void );
    PetscErrorCode ierr;

    // tmp1 = P * x
    Vec tmp1;
    ierr = VecDuplicate( x, &tmp1 );
    CHKERRQ( ierr );
    ierr = Apply_P( shellMat, x, tmp1 );
    CHKERRQ( ierr );

    // tmp2 = A * tmp1
    Vec tmp2;
    ierr = VecDuplicate( tmp1, &tmp2 );
    CHKERRQ( ierr );
    ierr = MatMult( ctx->A, tmp1, tmp2 );
    CHKERRQ( ierr );

    // y = P^T * tmp2  (P is symmetric in M-inner product; we implement simply P(tmp2))
    ierr = Apply_P( shellMat, tmp2, y );
    CHKERRQ( ierr );

    ierr = VecDestroy( &tmp1 );
    CHKERRQ( ierr );
    ierr = VecDestroy( &tmp2 );
    CHKERRQ( ierr );

    PetscFunctionReturn( 0 );
}

// y <- P^T * M * P * x
PetscErrorCode MP_mult( Mat shellMat, Vec x, Vec y )
{
    PetscFunctionBeginUser;
    void* ctx_void = nullptr;
    MatShellGetContext( shellMat, &ctx_void );
    ShellCtx* ctx = static_cast<ShellCtx*>( ctx_void );
    PetscErrorCode ierr;

    Vec tmp1, tmp2;
    ierr = VecDuplicate( x, &tmp1 );
    CHKERRQ( ierr );
    ierr = VecDuplicate( x, &tmp2 );
    CHKERRQ( ierr );

    ierr = Apply_P( shellMat, x, tmp1 );
    CHKERRQ( ierr ); // tmp1 = P x
    ierr = MatMult( ctx->M, tmp1, tmp2 );
    CHKERRQ( ierr ); // tmp2 = M * (P x)
    ierr = Apply_P( shellMat, tmp2, y );
    CHKERRQ( ierr ); // y = P^T tmp2 = P tmp2

    ierr = VecDestroy( &tmp1 );
    CHKERRQ( ierr );
    ierr = VecDestroy( &tmp2 );
    CHKERRQ( ierr );

    PetscFunctionReturn( 0 );
}

namespace assembly
{
int solveMixedEigenproblem( int argc,
                            char** argv,
                            const size_t nev,
                            const Eigen::SparseMatrix<double>& A_in,
                            const Eigen::SparseMatrix<double>& M_in,
                            const Eigen::SparseMatrix<double>& B_in )
{
    PetscErrorCode ierr;
    ierr = SlepcInitialize( &argc, &argv, nullptr, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );

    // Create matrices A, M, B (placeholders)
    Mat A = nullptr, M = nullptr, B = nullptr;
    // TODO: Decide local/global sizes for edge DOFs and node DOFs, and build mats accordingly.

    // For this skeleton we'll assume A, M are square (nedof x nedof) and B is (nodof x nedof).
    PetscInt nedof = A_in.cols(), nodof = B_in.rows(); // <-- fill these after mesh/FE assembly
    slepc_utils::eigenSparseToPetsc( A_in, A );
    slepc_utils::eigenSparseToPetsc( B_in, B );
    slepc_utils::eigenSparseToPetsc( M_in, M );

    PetscInt rows, cols;
    MatGetSize( A, &rows, &cols );
    std::cout << "A size: " << rows << " x " << cols << std::endl;
    MatGetSize( M, &rows, &cols );
    std::cout << "M size: " << rows << " x " << cols << std::endl;
    MatGetSize( B, &rows, &cols );
    std::cout << "B size: " << rows << " x " << cols << std::endl;
    // ---- ASSEMBLE A, M, B HERE ----
    // Example:
    // MatCreate(...); MatSetSizes(A, local_nedof, local_nedof, global_nedof, global_nedof); ...
    // Populate A, M, B with MatSetValues / MatSetValuesBlocked / MatAssemblyBegin/End
    //
    // After assembly set nedof, nodof accordingly.
    // --------------------------------


    // Create ShellCtx and inner KSPs
    ShellCtx ctx;
    ctx.A = A;
    ctx.M = M;
    ctx.B = B;
    if( !A || !M || !B )
    {
        PetscPrintf( PETSC_COMM_WORLD, "ERROR: You must assemble A, M, and B before running. Exiting.\n" );
        SlepcFinalize();
        return -1;
    }

    // Create Vecs sized appropriately (edge: nedof, nodal: nodof)
    Vec v_edge, v_edge_tmp;
    Vec v_node;
    ierr = MatCreateVecs( M, &v_edge, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // uses M's row layout
    // For tmp_Bv (nodal), create with B rows:
    ierr = MatCreateVecs( B, nullptr, &v_edge_tmp );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // b's col -> edge vec
    // Actually simpler: create node vec with B's rows
    ierr = MatCreateVecs( B, &v_node, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // node vec
    // But above MatCreateVecs usage depends on Mat shapes — adapt as necessary.

    // Create temporaries stored in ctx
    // Make them with correct communicators/sizes (use MatCreateVecs for M and B if helpful)
    // For clarity below we create all using MatCreateVecs with appropriate matrices:
    Vec tmp_edge, tmp_edge2;
    ierr = MatCreateVecs( M, &tmp_edge, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    ierr = MatCreateVecs( M, &tmp_edge2, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );

    ierr = MatCreateVecs( B, nullptr, &ctx.tmp_Bv );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // nodal
    ierr = MatCreateVecs( B, nullptr, &ctx.tmp_Sx );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // nodal
    ierr = MatCreateVecs( B, &ctx.tmp_Bt_x, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // edge (B has cols=edges)
    ierr = MatCreateVecs( M, &ctx.tmp_My2, nullptr );
    CHKERRABORT( PETSC_COMM_WORLD, ierr ); // edge

    // Create KSP for M solves (edge)
    KSP ksp_M;
    ierr = KSPCreate( PETSC_COMM_WORLD, &ksp_M );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    ierr = KSPSetOperators( ksp_M, M, M );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    // Name the KSP so you can tune via command line: -ksp_M_* options
    ierr = PetscObjectSetOptionsPrefix( (PetscObject)ksp_M, "ksp_M_" );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    KSPSetType( ksp_M, KSPCG );
    PC pc_M = nullptr;
    KSPGetPC( ksp_M, &pc_M );
    PCSetType( pc_M, PCJACOBI ); // or PCGAMG
    KSPSetFromOptions( ksp_M );

    // Create (or form) S = B * M^{-1} * B^T OR leave implicit and use ksp_S that applies S via MatMult
    Mat S = nullptr;
    // Create the shell matrix
    MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nodof, nodof, &ctx, &S);
    MatShellSetOperation(S, MATOP_MULT, (void(*)(void))Apply_S);

    // Store it in your context
    ctx.S_explicit = S;
    bool form_S_explicit = false; // set to true if you prefer to form S once (recommended if nodal DOFs small)
    if( form_S_explicit )
    {
        // Form S explicitly: S = B * (M^{-1} * B^T)
        // Strategy:
        // - For each column of B^T (i.e., each nodal basis vector e_i): compute y = B^T e_i (edge vec)
        // - Solve M x = y (use ksp_M as above)
        // - Then column i of S is B * x (nodal vec)
        // This is expensive but done once. For large problems prefer implicit.
        // Placeholder: user must implement this assembly. We create empty S with appropriate sizes.
        MatCreate( PETSC_COMM_WORLD, &S );
        MatSetSizes( S, PETSC_DECIDE, PETSC_DECIDE, nodof, nodof );
        MatSetFromOptions( S );
        MatSetUp( S );
        // Fill S entries here using the procedure above, then MatAssemblyBegin/End(S, MAT_FINAL_ASSEMBLY)
        ctx.S_is_explicit = true;
        ctx.S_explicit = S;
    }
    else
    {
        // If not forming S explicitly, you should implement a Mat that computes S*x := B * (M^{-1} * (B^T x))
        // For brevity we assume we will keep ksp_S solving against an explicit S; otherwise create another MatShell.
        ctx.S_is_explicit = false;
    }

    // Create KSP for S solves (nodal)
    KSP ksp_S;
    ierr = KSPCreate( PETSC_COMM_WORLD, &ksp_S );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    // if( ctx.S_is_explicit )
    {
        ierr = KSPSetOperators( ksp_S, ctx.S_explicit, ctx.S_explicit );
        CHKERRABORT( PETSC_COMM_WORLD, ierr );
    }
    // else
    // {
    //     // if S not explicit, user must create a MatShell for S and use it here
    //     PetscPrintf(
    //         PETSC_COMM_WORLD,
    //         "ERROR: implicit S operator not implemented in this skeleton. Form S explicitly or implement S-shell.\n" );
    //     SlepcFinalize();
    //     return -1;
    // }
    ierr = PetscObjectSetOptionsPrefix( (PetscObject)ksp_S, "ksp_S_" );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    KSPSetType( ksp_S, KSPCG );
    PC pc_S = nullptr;
    KSPGetPC( ksp_S, &pc_S );
    PCSetType( pc_S, PCNONE ); // or PCGAMG
    KSPSetFromOptions( ksp_S );

    // Store KSPs in ctx
    ctx.ksp_M = ksp_M;
    ctx.ksp_S = ksp_S;

    // Create MatShell objects for A_shell and M_shell
    Mat A_shell = nullptr, M_shell = nullptr;
    MatCreateShell( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nedof, nedof, &ctx, &A_shell );
    MatShellSetOperation( A_shell, MATOP_MULT, (void ( * )( void ))AP_mult );

    MatCreateShell( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nedof, nedof, &ctx, &M_shell );
    MatShellSetOperation( M_shell, MATOP_MULT, (void ( * )( void ))MP_mult );

    // Create EPS and set operators
    EPS eps;
    ierr = EPSCreate( PETSC_COMM_WORLD, &eps );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    ierr = EPSSetOperators( eps, A_shell, M_shell );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    ierr = EPSSetProblemType( eps, EPS_GHEP );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );

    // For shift-and-invert, we need to set target and which eigenpairs appropriately
    ierr = EPSSetWhichEigenpairs( eps, EPS_SMALLEST_REAL );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    
    // Set target for shift-and-invert (0 for largest magnitude)
    // ierr = EPSSetTarget( eps, 9.87 );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    
    // This is important for shift-and-invert with largest magnitude
    // ierr = EPSSetWhichEigenpairs( eps, EPS_TARGET_REAL );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );

    // Set number of eigenvalues to compute
    ierr = EPSSetDimensions( eps, nev, PETSC_DEFAULT, PETSC_DEFAULT );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    // Choose EPS type, tolerances, etc via command line or set here:
    // ierr = EPSSetType(eps, EPSLOBPCG); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    EPSSetFromOptions( eps );

    // Solve eigenproblem
    ierr = EPSSolve( eps );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );

    // Postprocess: get eigenpairs, verify residuals & constraints
    PetscInt nconv;
    ierr = EPSGetConverged( eps, &nconv );
    CHKERRABORT( PETSC_COMM_WORLD, ierr );
    PetscPrintf( PETSC_COMM_WORLD, "Number of converged eigenpairs: %d\n", (int)nconv );

    for( PetscInt i = 0; i < nconv; ++i )
    {
        Vec xr;
        VecCreate( PETSC_COMM_WORLD, &xr );
        VecSetFromOptions( xr );
        VecSetSizes( xr, PETSC_DECIDE, nedof );
        VecSetUp( xr );

        PetscScalar kr;
        ierr = EPSGetEigenpair( eps, i, &kr, nullptr, xr, nullptr ); // for generalized give xr, but M*xl may be needed
        CHKERRABORT( PETSC_COMM_WORLD, ierr );

        // Compute residual r = A*x - lambda*M*x
        Vec Ax, Mx;
        MatCreateVecs( A, &Ax, nullptr );
        MatCreateVecs( M, &Mx, nullptr );
        MatMult( A, xr, Ax );
        MatMult( M, xr, Mx );

        VecAXPY( Ax, -kr, Mx ); // Ax := Ax - kr*Mx => residual stored in Ax
        PetscReal resnorm;
        VecNorm( Ax, NORM_2, &resnorm );
        PetscPrintf( PETSC_COMM_WORLD,
                     "eig %d, lambda = %g, residual norm = %g\n",
                     (int)i,
                     (double)PetscRealPart( kr ),
                     (double)resnorm );

        // Compute constraint norm ||B*x||
        Vec Bx;
        MatCreateVecs( B, &Bx, nullptr ); // B rows -> nodal vec
        MatMult( B, xr, Bx );
        PetscReal bnorm;
        VecNorm( Bx, NORM_2, &bnorm );
        PetscPrintf( PETSC_COMM_WORLD, "  constraint ||B u|| = %g\n", (double)bnorm );

        // cleanup
        VecDestroy( &Bx );
        VecDestroy( &Ax );
        VecDestroy( &Mx );
        VecDestroy( &xr );
    }

    // cleanup: destroy objects (in practice destroy in reverse order)
    EPSDestroy( &eps );
    MatDestroy( &A_shell );
    MatDestroy( &M_shell );

    KSPDestroy( &ksp_M );
    KSPDestroy( &ksp_S );
    if( ctx.S_is_explicit ) MatDestroy( &ctx.S_explicit );

    // User-owned mats A,M,B should be destroyed by the user after they are done
    MatDestroy( &A );
    MatDestroy( &M );
    MatDestroy( &B );
    MatDestroy( &S );

    // destroy temporaries
    VecDestroy( &ctx.tmp_Bv );
    VecDestroy( &ctx.tmp_Sx );
    VecDestroy( &ctx.tmp_Bt_x );
    VecDestroy( &ctx.tmp_My2 );

    SlepcFinalize();
    return 0;
}
}