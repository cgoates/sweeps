#include <VTKOutput.hpp>
#include <Logging.hpp>
#include <NavierStokesDiscretization.hpp>
#include <SmallVector.hpp>
#include <CustomEigen.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <QuadratureRule.hpp>
#include <Eigen/LU>
#include <KnotVector.hpp>
#include <Eigen/Eigenvalues>

constexpr double SIGMA = 0;
constexpr double KINEMATIC_VISCOSITY = -1;

void localizeElement( api::NavierStokesDiscretization& nsd, const topology::Cell& elem )
{
    nsd.getH1().localizeElement( elem );
    nsd.getHDIV().localizeElement( elem );
    nsd.getL2().localizeElement( elem );
}

void localizePoint( api::NavierStokesDiscretization& nsd, const param::ParentPoint& ppt )
{
    nsd.getH1().localizePoint( ppt );
    nsd.getHDIV().localizePoint( ppt );
    nsd.getL2().localizePoint( ppt );
}

Eigen::MatrixXd localStiffness( api::NavierStokesDiscretization& nsd,
                                const assembly::QuadratureRule& quad_rule,
                                const size_t n_local_hdiv,
                                const size_t n_local_L2 )
{
    const double sigma = SIGMA;

    const double kinematic_viscosity = KINEMATIC_VISCOSITY;

    const size_t n_local_total = n_local_hdiv + n_local_L2;

    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero( n_local_total, n_local_total );

    const Eigen::PermutationMatrix<4> symmetric_permutation( Eigen::Vector4i( 0, 2, 1, 3 ));

    quad_rule.iterateQuadraturePoints( [&]( const param::ParentPoint& qpt, const double qwt ) {
        localizePoint( nsd, qpt );
        const double jac_det = nsd.getH1().evaluateJacobian( nsd.controlPoints() ).determinant();

        const Eigen::MatrixXd transformed_basis =
            eval::piolaTransformedHDivBasis( nsd.getHDIV(), nsd.getH1(), nsd.controlPoints() );
        const Eigen::VectorXd transformed_basis_L2 =
            eval::piolaTransformedL2Basis( nsd.getL2(), nsd.getH1(), nsd.controlPoints() );
        const Eigen::MatrixXd gradient_basis =
            eval::piolaTransformedHDivFirstDerivatives( nsd.getHDIV(), nsd.getH1(), nsd.controlPoints() );

        const Eigen::MatrixXd sym_gradients = 0.5 * ( gradient_basis + gradient_basis * symmetric_permutation );

        // Viscous and pressure contributions
        for( size_t a = 0; a < n_local_hdiv; ++a )
        {
            const auto grad_a = sym_gradients.row(a);//FIXME: This is gradient in parametric space
            const auto basis_a = transformed_basis.row(a);
            for( size_t b = 0; b < n_local_hdiv; ++b )
            {
                const auto grad_b = sym_gradients.row(b);//FIXME: This is gradient in parametric space
                const auto basis_b = transformed_basis.row(b);

                const double term1 = 2 * kinematic_viscosity * grad_a.dot( grad_b );
                const double term2 = sigma * basis_a.dot( basis_b );

                ke( a, b ) += ( term1 + term2 ) * jac_det * qwt;
            }
        }

        // Pressure coupling
        for( size_t a = 0; a < n_local_hdiv; ++a )
        {
            const double trace_grad = gradient_basis( a, 0 ) + gradient_basis( a, 3 );
            for( size_t b = 0; b < n_local_L2; ++b )
            {
                const double ke_val = -trace_grad * transformed_basis_L2( b ) * jac_det * qwt;
                ke( a, b + n_local_hdiv ) += ke_val;
                ke( b + n_local_hdiv, a ) += -ke_val;
            }
        }
    } );

    return ke;
}

Eigen::SparseMatrix<double> assembleStiffnessMatrix( api::NavierStokesDiscretization& nsd, const assembly::QuadratureRule& quad_rule )
{
    // # global SIGMA  # Use global SIGMA
    // # SIGMA = 1  # Set for Stokes flow

    const size_t n_hdiv = nsd.getHDIV().splineSpace().numFunctions();
    const size_t n_l2 = nsd.getL2().splineSpace().numFunctions();
    const size_t n_total_funcs = n_hdiv + n_l2;
    const topology::CombinatorialMap& cmap = nsd.getHDIV().splineSpace().basisComplex().parametricAtlas().cmap();

    std::vector<Eigen::Triplet<double>> triplets;

    Eigen::VectorXi ID = Eigen::VectorXi::LinSpaced( n_total_funcs, 0, n_total_funcs - 1 );

    // ID( n_hdiv ) = -1;  // Mark first pressure term as known (Dirichlet boundary)

    topology::iterateCellsWhile( cmap, cmap.dim(), [&]( const topology::Cell& elem ) {
        localizeElement( nsd, elem );
        const std::vector<basis::FunctionId> elem_connectivity_HDIV = nsd.getHDIV().splineSpace().connectivity( elem );
        const std::vector<basis::FunctionId> elem_connectivity_L2 = nsd.getL2().splineSpace().connectivity( elem );
        const size_t n_local_hdiv = elem_connectivity_HDIV.size();
        const size_t n_local_L2 = elem_connectivity_L2.size();

        // Use Stokes stiffness since SIGMA=1
        const Eigen::MatrixXd k_e = localStiffness( nsd, quad_rule, n_local_hdiv, n_local_L2 );

        for( size_t a = 0; a < n_local_hdiv; ++a )
        {
            const Eigen::Index A = elem_connectivity_HDIV.at( a ).id();
            if( ID( A ) != -1 )
            {
                // F( A ) += f_e( a );
                for( size_t b = 0; b < n_local_hdiv; ++b )
                {
                    const Eigen::Index B = elem_connectivity_HDIV.at( b ).id();
                    if( ID( B ) != -1 )
                    {
                        triplets.emplace_back( A, B, k_e( a, b ) );
                    }
                }
                for( size_t b = 0; b < n_local_L2; ++b )
                {
                    const Eigen::Index B = elem_connectivity_L2.at( b ).id() + n_hdiv;
                    if( ID( B ) != -1 )
                    {
                        triplets.emplace_back( A, B, k_e( a, b + n_local_hdiv ) );
                        triplets.emplace_back( B, A, k_e( b + n_local_hdiv, a ) );
                    }
                }
            }
        }
        return true;
    } );

    Eigen::SparseMatrix<double> K( n_total_funcs, n_total_funcs );
    K.setFromTriplets(triplets.begin(), triplets.end());

    return K;
}

Eigen::SparseMatrix<double> pressureMassMatrix( api::NavierStokesDiscretization& nsd, const assembly::QuadratureRule& quad_rule )
{
    const size_t n_l2 = nsd.getL2().splineSpace().numFunctions();
    const topology::CombinatorialMap& cmap = nsd.getL2().splineSpace().basisComplex().parametricAtlas().cmap();

    std::vector<Eigen::Triplet<double>> triplets;

    topology::iterateCellsWhile( cmap, cmap.dim(), [&]( const topology::Cell& elem ) {
        localizeElement( nsd, elem );
        const std::vector<basis::FunctionId> elem_connectivity_L2 = nsd.getL2().splineSpace().connectivity( elem );
        const size_t n_local_L2 = elem_connectivity_L2.size();

        Eigen::MatrixXd me = Eigen::MatrixXd::Zero( n_local_L2, n_local_L2 );

        quad_rule.iterateQuadraturePoints( [&]( const param::ParentPoint& qpt, const double qwt ) {
            localizePoint( nsd, qpt );
            const double jac_det = nsd.getH1().evaluateJacobian( nsd.controlPoints() ).determinant();

            const Eigen::VectorXd transformed_basis_L2 =
                eval::piolaTransformedL2Basis( nsd.getL2(), nsd.getH1(), nsd.controlPoints() );
            
            for( size_t a = 0; a < n_local_L2; ++a )
            {
                const double basis_a = transformed_basis_L2( a );
                for( size_t b = 0; b < n_local_L2; ++b )
                {
                    const double basis_b = transformed_basis_L2( b );
                    me( a, b ) += basis_a * basis_b * jac_det * qwt;
                }
            }
        } );

        for( size_t a = 0; a < n_local_L2; ++a )
        {
            const Eigen::Index A = elem_connectivity_L2.at( a ).id();
            for( size_t b = 0; b < n_local_L2; ++b )
            {
                const Eigen::Index B = elem_connectivity_L2.at( b ).id();
                triplets.emplace_back( A, B, me( a, b ) );
            }
        }
        return true;
    } );

    Eigen::SparseMatrix<double> Q( n_l2, n_l2 );
    Q.setFromTriplets(triplets.begin(), triplets.end());
    return Q;
}

Eigen::SparseMatrix<double> matrixN( api::NavierStokesDiscretization& nsd,
                                     const assembly::QuadratureRule& quad_rule,
                                     const Eigen::SparseMatrix<double>& K )
{
    const Eigen::Index n_hdiv = nsd.getHDIV().splineSpace().numFunctions();
    const Eigen::SparseMatrix<double> Q = pressureMassMatrix( nsd, quad_rule );
    std::cout << "Q( " << Q.rows() << ", " << Q.cols() << " ):" << std::endl;
    std::cout << Q.toDense() << std::endl;

    Eigen::SparseMatrix<double> N( K.rows(), K.cols() );
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve( K.nonZeros() );
    for( Eigen::Index i = 0; i < n_hdiv; i++ )
    {
        for( Eigen::SparseMatrix<double>::InnerIterator it( K, i ); it and it.row() < n_hdiv; ++it )
        {
            triplets.emplace_back( it.row(), it.col(), it.value() );
        }
    }

    for( Eigen::Index i = 0; i < Q.cols(); i++ )
    {
        for( Eigen::SparseMatrix<double>::InnerIterator it( Q, i ); it; ++it )
        {
            triplets.emplace_back( it.row() + n_hdiv, it.col() + n_hdiv, it.value() );
        }
    }
    N.setFromTriplets( triplets.begin(), triplets.end() );
    return N;
}

int main()
{
    const size_t degree = 2;
    const basis::KnotVector kv_s( { 0.0, 0.0, 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.0, 1.0 }, 1e-9 );
    const basis::KnotVector kv_t = kv_s;
    const param::ParentDomain domain = param::cubeDomain( 2 );
    const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( degree + 1, domain );
    api::NavierStokesTPDiscretization nsd(
        kv_s,
        kv_t,
        degree,
        degree,
        util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ) } ).transpose() );

    const Eigen::SparseMatrix<double> K = assembleStiffnessMatrix( nsd, quad_rule );
    std::cout << "K( " << K.rows() << ", " << K.cols() << " ):" << std::endl;
    std::cout << K.toDense() << std::endl;
    const Eigen::SparseMatrix<double> N = matrixN( nsd, quad_rule, K );
    std::cout << "N( " << N.rows() << ", " << N.cols() << " ):" << std::endl;
    std::cout << N.toDense() << std::endl;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver( K.toDense(), N.toDense() );
    std::cout << "Eigenvalues:" << std::endl;
    std::cout << solver.eigenvalues() << std::endl;
    std::cout << "Eigenvectors (by column):" << std::endl;
    std::cout << solver.eigenvectors() << std::endl;
    return 0;
}
