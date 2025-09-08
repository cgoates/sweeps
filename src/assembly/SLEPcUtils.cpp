#include <slepceps.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <tuple>

namespace slepc_utils
{
    void checkErrorCode( PetscErrorCode ierr )
    {
        if( ierr ) throw std::runtime_error( "PETSc error code: " + std::to_string( ierr ) );
    }

    // Convert Eigen sparse matrix to PETSc matrix
    PetscErrorCode eigenSparseToPetsc( const Eigen::SparseMatrix<double>& eigenMat, Mat& petscMat, MPI_Comm comm )
    {
        PetscErrorCode ierr;
        PetscInt m = eigenMat.rows();
        PetscInt n = eigenMat.cols();

        ierr = MatCreate( comm, &petscMat );
        checkErrorCode( ierr );
        ierr = MatSetSizes( petscMat, PETSC_DECIDE, PETSC_DECIDE, m, n );
        checkErrorCode( ierr );
        ierr = MatSetFromOptions( petscMat );
        checkErrorCode( ierr );

        // Get local ownership range
        PetscInt rstart, rend;
        ierr = MatGetOwnershipRange( petscMat, &rstart, &rend );
        checkErrorCode( ierr );

        // Set values from sparse matrix
        for( PetscInt k = 0; k < eigenMat.cols(); k++ )
        {
            for( Eigen::SparseMatrix<double>::InnerIterator it( eigenMat, k ); it; ++it )
            {
                PetscInt col = it.col();
                PetscInt row = it.row();
                if( row < rstart || row >= rend ) continue; // Only set values for owned rows
                PetscScalar val = it.value();
                ierr = MatSetValue( petscMat, row, col, val, INSERT_VALUES );
                checkErrorCode( ierr );
            }
        }

        ierr = MatAssemblyBegin( petscMat, MAT_FINAL_ASSEMBLY );
        checkErrorCode( ierr );
        ierr = MatAssemblyEnd( petscMat, MAT_FINAL_ASSEMBLY );
        checkErrorCode( ierr );

        return 0;
    }

    // Convert PETSc vector to Eigen vector
    Eigen::VectorXd petscToEigen( const Vec& petscVec, MPI_Comm comm = PETSC_COMM_WORLD )
    {
        PetscErrorCode ierr;
        PetscInt n;

        ierr = VecGetSize( petscVec, &n );
        checkErrorCode( ierr );
        Eigen::VectorXd eigenVec( n );

        // Get local ownership range
        PetscInt rstart, rend;
        ierr = VecGetOwnershipRange( petscVec, &rstart, &rend );
        checkErrorCode( ierr );

        // Get array of values
        const PetscScalar* array;
        ierr = VecGetArrayRead( petscVec, &array );
        checkErrorCode( ierr );

        for( PetscInt i = rstart; i < rend; i++ )
        {
            eigenVec( i ) = array[i - rstart];
        }

        ierr = VecRestoreArrayRead( petscVec, &array );
        checkErrorCode( ierr );

        // Gather all values to all processes
        MPI_Allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, eigenVec.data(), NULL, NULL, MPIU_SCALAR, comm );

        return eigenVec;
    }

PetscErrorCode createK1Matrix( const Mat& Q, const Mat& K1_bar, Mat& A, MPI_Comm comm = PETSC_COMM_WORLD )
{
    PetscErrorCode ierr;

    // Create KSP solver for Q
    KSP ksp;
    ierr = KSPCreate( comm, &ksp );
    checkErrorCode( ierr );
    ierr = KSPSetOperators( ksp, Q, Q );
    checkErrorCode( ierr );

    // Use iterative solver with preconditioner
    PC pc = nullptr;
    ierr = KSPGetPC( ksp, &pc );
    checkErrorCode( ierr );
    ierr = PCSetType( pc, PCGAMG );
    checkErrorCode( ierr );
    ierr = KSPSetType( ksp, KSPGMRES );
    checkErrorCode( ierr );
    ierr = KSPSetTolerances( ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    checkErrorCode( ierr );

    ierr = KSPSetUp( ksp );
    checkErrorCode( ierr );

    // Get matrix dimensions
    PetscInt m, n_rhs;
    ierr = MatGetSize(K1_bar, &m, &n_rhs);
    checkErrorCode( ierr );

    // Create solution matrix
    Mat X;
    ierr = MatDuplicate( K1_bar, MAT_DO_NOT_COPY_VALUES, &X );
    checkErrorCode( ierr );

    // Create vectors for individual RHS and solution
    Vec rhs, sol;
    ierr = MatCreateVecs(Q, &rhs, &sol);
    checkErrorCode( ierr );

    // Get ownership range for parallel processing
    PetscInt Istart, Iend;
    ierr = MatGetOwnershipRange(K1_bar, &Istart, &Iend);
    checkErrorCode( ierr );

    // Solve for each column individually
    for (PetscInt col = 0; col < n_rhs; col++) {
        // Get column from K1_bar as RHS vector
        ierr = MatGetColumnVector(K1_bar, rhs, col);
        checkErrorCode( ierr );
        
        // Solve Q * x = rhs
        ierr = KSPSolve(ksp, rhs, sol);
        checkErrorCode( ierr );
        
        // Get local portion of solution vector
        const PetscScalar *sol_array;
        ierr = VecGetArrayRead(sol, &sol_array);
        checkErrorCode( ierr );
        
        // Set values in solution matrix X
        PetscInt nlocal;
        ierr = VecGetLocalSize(sol, &nlocal);
        checkErrorCode( ierr );
        
        PetscInt *indices = new PetscInt[nlocal];
        for (PetscInt i = 0; i < nlocal; i++) {
            indices[i] = Istart + i;
        }
        
        ierr = MatSetValues(X, nlocal, indices, 1, &col, sol_array, INSERT_VALUES);
        checkErrorCode( ierr );
        
        delete[] indices;
        ierr = VecRestoreArrayRead(sol, &sol_array);
        checkErrorCode( ierr );
    }

    // Assemble the solution matrix
    ierr = MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY);
    checkErrorCode( ierr );
    ierr = MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY);
    checkErrorCode( ierr );

    // Now compute A = K1_bar^T * X
    ierr = MatTransposeMatMult( K1_bar, X, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A );
    checkErrorCode( ierr );

    // Cleanup
    ierr = MatDestroy( &X );
    checkErrorCode( ierr );
    ierr = VecDestroy( &rhs );
    checkErrorCode( ierr );
    ierr = VecDestroy( &sol );
    checkErrorCode( ierr );
    ierr = KSPDestroy( &ksp );
    checkErrorCode( ierr );

    return 0;
}

std::pair<std::vector<PetscScalar>, std::vector<Eigen::VectorXd>>
    solveGeneralizedEigenvalueSparse( const Eigen::SparseMatrix<double>& K1_bar,
                                      const Eigen::SparseMatrix<double>& Q,
                                      const Eigen::SparseMatrix<double>& K2,
                                      int nev,
                                      MPI_Comm comm )
{
    PetscErrorCode ierr;

    // Convert matrices to PETSc format
    Mat Q_petsc, K1_bar_petsc, K2_petsc;
    ierr = eigenSparseToPetsc( Q, Q_petsc, comm );
    checkErrorCode( ierr );
    ierr = eigenSparseToPetsc( K1_bar, K1_bar_petsc, comm );
    checkErrorCode( ierr );
    ierr = eigenSparseToPetsc( K2, K2_petsc, comm );
    checkErrorCode( ierr );

    // std::cout << K1_bar << std::endl << std::endl;
    // MatView( K1_bar_petsc, PETSC_VIEWER_STDOUT_WORLD );

    // std::cout << "\n---------\n" << Q << std::endl << std::endl;
    // MatView( Q_petsc, PETSC_VIEWER_STDOUT_WORLD );
    // std::cout  << "\n---------\n" << K2 << std::endl << std::endl;
    // MatView( K2_petsc, PETSC_VIEWER_STDOUT_WORLD );

    Mat K1;
    ierr = createK1Matrix( Q_petsc, K1_bar_petsc, K1, comm );
    checkErrorCode( ierr );

    // Create EPS solver
    EPS eps;
    ierr = EPSCreate( comm, &eps );
    checkErrorCode( ierr );

    // Set the problem type
    ierr = EPSSetOperators( eps, K2_petsc, K1 );
    checkErrorCode( ierr );
    // ierr = EPSSetProblemType( eps, EPS_GHEP );//This will be automatically detected if not set
    // checkErrorCode( ierr );

    // For shift-and-invert, we need to set target and which eigenpairs appropriately
    ierr = EPSSetWhichEigenpairs( eps, EPS_SMALLEST_MAGNITUDE );
    checkErrorCode( ierr );
    
    // Set target for shift-and-invert (0 for largest magnitude)
    ierr = EPSSetTarget( eps, 0.0 );
    checkErrorCode( ierr );
    
    // This is important for shift-and-invert with largest magnitude
    ierr = EPSSetWhichEigenpairs( eps, EPS_TARGET_MAGNITUDE );
    checkErrorCode( ierr );

    // Set number of eigenvalues to compute
    ierr = EPSSetDimensions( eps, nev, PETSC_DEFAULT, PETSC_DEFAULT );
    checkErrorCode( ierr );

    // Configure iterative solver for shift-and-invert
    ST st;
    KSP ksp;
    PC pc;
    
    ierr = EPSGetST(eps, &st);
    checkErrorCode(ierr);
    ierr = STSetType(st, STSINVERT);
    checkErrorCode(ierr);
    ierr = STGetKSP(st, &ksp);
    checkErrorCode(ierr);
    ierr = KSPSetType(ksp, KSPCG);
    checkErrorCode(ierr);
    ierr = KSPGetPC(ksp, &pc);
    checkErrorCode(ierr);
    ierr = PCSetType(pc, PCJACOBI);
    checkErrorCode(ierr);
    ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    checkErrorCode(ierr);

    // Set solver options from command line
    ierr = EPSSetFromOptions( eps );
    checkErrorCode( ierr );

    // Solve the eigenvalue problem
    ierr = EPSSolve( eps );
    checkErrorCode( ierr );

    // Get number of converged eigenvalues
    PetscInt nconv;
    ierr = EPSGetConverged( eps, &nconv );
    checkErrorCode( ierr );

    std::cout << "Number of converged eigenvalues: " << nconv << std::endl;

    // Extract eigenvalues and eigenvectors
    std::vector<PetscScalar> eigenvalues( nconv );
    std::vector<Eigen::VectorXd> eigenvectors( nconv );

    for( PetscInt i = 0; i < nconv; i++ )
    {
        Vec eigenvector;
        ierr = MatCreateVecs( K2_petsc, &eigenvector, NULL );
        checkErrorCode( ierr );

        PetscScalar eigenvalue;
        ierr = EPSGetEigenpair( eps, i, &eigenvalue, NULL, eigenvector, NULL );
        checkErrorCode( ierr );

        eigenvalues[i] = eigenvalue;
        eigenvectors[i] = petscToEigen( eigenvector, comm );

        ierr = VecDestroy( &eigenvector );
        checkErrorCode( ierr );
    }

    // Clean up
    ierr = EPSDestroy( &eps );
    checkErrorCode( ierr );
    ierr = MatDestroy( &Q_petsc );
    checkErrorCode( ierr );
    ierr = MatDestroy( &K1_bar_petsc );
    checkErrorCode( ierr );
    ierr = MatDestroy( &K2_petsc );
    checkErrorCode( ierr );
    ierr = MatDestroy( &K1 );
    checkErrorCode( ierr );

    return std::make_pair( eigenvalues, eigenvectors );
}

std::pair<std::vector<PetscScalar>, std::vector<Eigen::VectorXd>>
    solveGeneralizedEigenvalueSparse( const Eigen::SparseMatrix<double>& A,
                                      const Eigen::SparseMatrix<double>& B,
                                      int nev,
                                      int n_edge,
                                      MPI_Comm comm )
{
    PetscErrorCode ierr;

    const int n_node = A.rows() - n_edge;

    // Convert matrices to PETSc format
    Mat B_petsc, A_petsc;
    ierr = eigenSparseToPetsc( B, B_petsc, comm );
    checkErrorCode( ierr );
    ierr = eigenSparseToPetsc( A, A_petsc, comm );
    checkErrorCode( ierr );

    // Create EPS solver
    EPS eps;
    ierr = EPSCreate( comm, &eps );
    checkErrorCode( ierr );

    // Set the problem type
    ierr = EPSSetOperators( eps, A_petsc, B_petsc );
    checkErrorCode( ierr );
    ierr = EPSSetProblemType( eps, EPS_GHEP );//This will be automatically detected if not set
    // checkErrorCode( ierr );

    // For shift-and-invert, we need to set target and which eigenpairs appropriately
    ierr = EPSSetWhichEigenpairs( eps, EPS_LARGEST_MAGNITUDE );
    checkErrorCode( ierr );
    
    // Set target for shift-and-invert (0 for largest magnitude)
    // ierr = EPSSetTarget( eps, 0.0 );
    checkErrorCode( ierr );
    
    // This is important for shift-and-invert with largest magnitude
    ierr = EPSSetWhichEigenpairs( eps, EPS_TARGET_MAGNITUDE );
    checkErrorCode( ierr );

    // Set number of eigenvalues to compute
    ierr = EPSSetDimensions( eps, nev, PETSC_DEFAULT, PETSC_DEFAULT );
    checkErrorCode( ierr );

    // Configure iterative solver for shift-and-invert
    ST st;
    KSP ksp;
    PC pc;

    IS is_edge, is_node;
    ISCreateStride(PETSC_COMM_WORLD, n_edge, 0, 1, &is_edge);
    ISCreateStride(PETSC_COMM_WORLD, n_node, n_edge, 1, &is_node);
    
    ierr = EPSGetST(eps, &st);
    checkErrorCode(ierr);
    ierr = STSetType(st, STSHIFT);
    ierr = STSetShift( st, 1e-6 );
    checkErrorCode(ierr);
    ierr = STGetKSP(st, &ksp);
    checkErrorCode(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY);
    checkErrorCode(ierr);
    ierr = KSPGetPC(ksp, &pc);
    checkErrorCode(ierr);
    ierr = PCSetType(pc, PCFIELDSPLIT);
    ierr = PCFieldSplitSetIS(pc, "edge", is_edge);
    ierr = PCFieldSplitSetIS(pc, "node", is_node);
    ierr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
    checkErrorCode(ierr);
    ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    checkErrorCode(ierr);

    // Set solver options from command line
    ierr = EPSSetFromOptions( eps );
    checkErrorCode( ierr );

    // Solve the eigenvalue problem
    ierr = EPSSolve( eps );
    checkErrorCode( ierr );

    // Get number of converged eigenvalues
    PetscInt nconv;
    ierr = EPSGetConverged( eps, &nconv );
    checkErrorCode( ierr );

    std::cout << "Number of converged eigenvalues: " << nconv << std::endl;

    // Extract eigenvalues and eigenvectors
    std::vector<PetscScalar> eigenvalues( nconv );
    std::vector<Eigen::VectorXd> eigenvectors( nconv );

    for( PetscInt i = 0; i < nconv; i++ )
    {
        Vec eigenvector;
        ierr = MatCreateVecs( A_petsc, &eigenvector, NULL );
        checkErrorCode( ierr );

        PetscScalar eigenvalue;
        PetscScalar eigenvalue_im;
        ierr = EPSGetEigenpair( eps, i, &eigenvalue, &eigenvalue_im, eigenvector, NULL );
        checkErrorCode( ierr );

        std::cout << "( " << eigenvalue << ", " << eigenvalue_im << " ), ";

        eigenvalues[i] = eigenvalue;
        eigenvectors[i] = petscToEigen( eigenvector, comm );

        ierr = VecDestroy( &eigenvector );
        checkErrorCode( ierr );
    }
    std::cout << std::endl;
    ISDestroy( &is_edge );
    ISDestroy( &is_node );

    // Clean up
    ierr = EPSDestroy( &eps );
    checkErrorCode( ierr );
    ierr = MatDestroy( &B_petsc );
    checkErrorCode( ierr );
    ierr = MatDestroy( &A_petsc );
    checkErrorCode( ierr );

    return std::make_pair( eigenvalues, eigenvectors );
}

} // namespace slepc_utils