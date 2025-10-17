#include <VTKOutput.hpp>
#include <Logging.hpp>
#include <DeRhamComplexDiscretization.hpp>
#include <SmallVector.hpp>
#include <CustomEigen.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <QuadratureRule.hpp>
#include <Eigen/LU>
#include <KnotVector.hpp>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <numbers>
#include <CommonUtils.hpp>
#include <NavierStokesDiscretization.hpp>
#include <SLEPcUtils.hpp>
#include <SaddlePointSolve.hpp>
#include <unsupported/Eigen/SparseExtra>
#include <DeRhamComplexTestCases.hpp>

using namespace assembly;

Eigen::MatrixXd localStiffness( DeRhamComplexDiscretization& drcd,
                                const assembly::QuadratureRule& quad_rule,
                                const size_t n_local_hcurl,
                                const size_t n_local_h1 )
{
    const size_t n_local_total = n_local_hcurl + n_local_h1;
    const size_t spatial_dim = drcd.controlPoints().rows();

    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero( n_local_total, n_local_total );

    quad_rule.iterateQuadraturePoints( [&]( const param::ParentPoint& qpt, const double qwt ) {
        localizePoint( drcd, qpt );
        const double jac_det = drcd.getManifold().evaluateJacobian( drcd.controlPoints() ).determinant();

        const Eigen::MatrixXd transformed_basis_hcurl =
            // drcd.getHCURL().evaluateBasis();
            eval::piolaTransformedHCurlBasis( drcd.getHCURL(), drcd.getManifold(), drcd.controlPoints() );
        const Eigen::MatrixXd derivative_basis_hcurl =
            // drcd.getHCURL().evaluateFirstDerivativesFromParamToSpatial();
            eval::piolaTransformedHCurlFirstDerivatives( drcd.getHCURL(), drcd.getManifold(), drcd.controlPoints() );
        const Eigen::MatrixXd derivative_basis_h1 = 
            // drcd.getH1().evaluateFirstDerivativesFromParamToSpatial();
            eval::piolaTransformedH1FirstDerivatives( drcd.getH1(), drcd.getManifold(), drcd.controlPoints() );

        const Eigen::MatrixXd curl_hcurl = [&]() -> Eigen::MatrixXd {
            if( spatial_dim == 2 )
            {
                return derivative_basis_hcurl.col( 1 ) - derivative_basis_hcurl.col( 2 );
            }
            else if( spatial_dim == 3 )
            {
                return ( Eigen::MatrixXd( n_local_hcurl, 3 ) <<
                    derivative_basis_hcurl.col( 5 ) - derivative_basis_hcurl.col( 7 ), // dNz/dy - dNy/dz
                    derivative_basis_hcurl.col( 6 ) - derivative_basis_hcurl.col( 2 ), // dNx/dz - dNz/dx
                    derivative_basis_hcurl.col( 1 ) - derivative_basis_hcurl.col( 3 ) // dv/dx - du/dy
                ).finished();
            }
            else
                throw std::runtime_error( "Unsupported spatial dimension for curl computation." );
        }();

        // HCURL x HCURL block
        for( size_t a = 0; a < n_local_hcurl; ++a )
        {
            const auto curl_a = curl_hcurl.row( a );
            for( size_t b = 0; b < n_local_hcurl; ++b )
            {
                const auto curl_b = curl_hcurl.row( b );

                ke( a, b ) += curl_a.dot( curl_b ) * jac_det * qwt;
            }
        }

        // std::cout << "transformed basis HCURL at qpt\n" << transformed_basis_hcurl << std::endl;
        // std::cout << "derivative basis H1 at qpt\n" << derivative_basis_h1 << std::endl;

        // HCURL x H1 and H1 x HCURL blocks
        for( size_t a = 0; a < n_local_hcurl; ++a )
        {
            const auto val_a = transformed_basis_hcurl.row( a );
            for( size_t b = 0; b < n_local_h1; ++b )
            {
                const double ke_val = val_a.dot( derivative_basis_h1.row( b ) ) * jac_det * qwt;
                ke( a, b + n_local_hcurl ) += ke_val;
                ke( b + n_local_hcurl, a ) += ke_val;
            }
        }
    } );

    return ke;
}

Eigen::SparseMatrix<double> assembleStiffnessMatrix( DeRhamComplexDiscretization& drcd, const assembly::QuadratureRule& quad_rule, const Eigen::VectorXi& id_array )
{
    const size_t n_constraints = ( id_array.array() == -1 ).count();
    const size_t n_hcurl = drcd.getHCURL().splineSpace().numFunctions();
    const size_t n_h1 = drcd.getH1().splineSpace().numFunctions();
    const size_t n_total_funcs = n_hcurl + n_h1;
    const topology::CombinatorialMap& cmap = drcd.getHCURL().splineSpace().basisComplex().parametricAtlas().cmap();

    std::vector<Eigen::Triplet<double>> triplets;

    topology::iterateCellsWhile( cmap, cmap.dim(), [&]( const topology::Cell& elem ) {
        localizeElement( drcd, elem );
        const std::vector<basis::FunctionId> elem_connectivity_HCURL = drcd.getHCURL().splineSpace().connectivity( elem );
        const std::vector<basis::FunctionId> elem_connectivity_H1 = drcd.getH1().splineSpace().connectivity( elem );
        const size_t n_local_hcurl = elem_connectivity_HCURL.size();
        const size_t n_local_H1 = elem_connectivity_H1.size();

        const Eigen::MatrixXd k_e = localStiffness( drcd, quad_rule, n_local_hcurl, n_local_H1 );
        // std::cout << "element " << elem << " stiffness matrix:\n" << k_e << std::endl;

        for( size_t a = 0; a < n_local_hcurl; ++a )
        {
            const Eigen::Index A = elem_connectivity_HCURL.at( a ).id();
            if( id_array( A ) == -1 )
                continue;
            for( size_t b = 0; b < n_local_hcurl; ++b )
            {
                const Eigen::Index B = elem_connectivity_HCURL.at( b ).id();
                if( id_array( B ) == -1 )
                    continue;
                triplets.emplace_back( id_array( A ), id_array( B ), k_e( a, b ) );
            }

            for( size_t b = 0; b < n_local_H1; ++b )
            {
                const Eigen::Index B = elem_connectivity_H1.at( b ).id() + n_hcurl;
                if( id_array( B ) == -1 )
                    continue;
                // std::cout << "a,b: " << a << ", " << b << " B: " << B << " " << k_e( a, b + n_local_hcurl ) << std::endl;
                triplets.emplace_back( id_array( A ), id_array( B ), k_e( a, b + n_local_hcurl ) );
                triplets.emplace_back( id_array( B ), id_array( A ), k_e( b + n_local_hcurl, a ) );
            }
        }
        return true;
    } );

    Eigen::SparseMatrix<double> K( n_total_funcs - n_constraints, n_total_funcs - n_constraints );
    K.setFromTriplets(triplets.begin(), triplets.end());

    return K;
}

Eigen::SparseMatrix<double> massMatrix( DeRhamComplexDiscretization& drcd, const assembly::QuadratureRule& quad_rule, const Eigen::VectorXi& id_array )
{
    const size_t n_constraints = ( id_array.array() == -1 ).count();
    const size_t n_total_funcs = drcd.getHCURL().splineSpace().numFunctions() + drcd.getH1().splineSpace().numFunctions();
    std::cout << "Total number of functions: " << n_total_funcs << " of which " << drcd.getHCURL().splineSpace().numFunctions() << " are in HCURL." << std::endl;
    const topology::CombinatorialMap& cmap = drcd.getHCURL().splineSpace().basisComplex().parametricAtlas().cmap();

    std::vector<Eigen::Triplet<double>> triplets;

    topology::iterateCellsWhile( cmap, cmap.dim(), [&]( const topology::Cell& elem ) {
        localizeElement( drcd, elem );
        const std::vector<basis::FunctionId> elem_connectivity_HCURL = drcd.getHCURL().splineSpace().connectivity( elem );
        const size_t n_local_HCURL = elem_connectivity_HCURL.size();

        Eigen::MatrixXd me = Eigen::MatrixXd::Zero( n_local_HCURL, n_local_HCURL );

        quad_rule.iterateQuadraturePoints( [&]( const param::ParentPoint& qpt, const double qwt ) {
            localizePoint( drcd, qpt );
            const double jac_det = drcd.getManifold().evaluateJacobian( drcd.controlPoints() ).determinant();

            const Eigen::MatrixXd transformed_basis_HCURL =
                // drcd.getHCURL().evaluateBasis();
                eval::piolaTransformedHCurlBasis( drcd.getHCURL(), drcd.getManifold(), drcd.controlPoints() );
            
            for( size_t a = 0; a < n_local_HCURL; ++a )
            {
                const auto basis_a = transformed_basis_HCURL.row( a );
                for( size_t b = 0; b < n_local_HCURL; ++b )
                {
                    const auto basis_b = transformed_basis_HCURL.row( b );
                    me( a, b ) += basis_a.dot( basis_b ) * jac_det * qwt;
                }
            }
        } );
        // std::cout << "n_local_HCURL: " << n_local_HCURL << std::endl;
        // std::cout << "Element mass matrix:\n" << me << std::endl;

        for( size_t a = 0; a < n_local_HCURL; ++a )
        {
            const Eigen::Index A = elem_connectivity_HCURL.at( a ).id();
            if( id_array( A ) == -1 )
                continue;
            for( size_t b = 0; b < n_local_HCURL; ++b )
            {
                const Eigen::Index B = elem_connectivity_HCURL.at( b ).id();
                if( id_array( B ) == -1 )
                    continue;
                triplets.emplace_back( id_array( A ), id_array( B ), me( a, b ) );
            }
        }
        return true;
    } );

    Eigen::SparseMatrix<double> M( n_total_funcs - n_constraints, n_total_funcs - n_constraints );
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

Eigen::VectorXi idArray( const DeRhamComplexTPDiscretization& drcd )
{
    const auto& hcurl = drcd.HCURL_ss;
    const auto& component_bases = hcurl.scalarTPBases();
    const auto tangent_on_side = [&]( const api::PatchSide& side ) {
        const bool is_S = side == api::PatchSide::S0 or side == api::PatchSide::S1;

        const util::IndexVec lengths = getTPLengths( is_S ? *component_bases.at( 1 ) : *component_bases.at( 0 ) );

        const auto iter_dir =
            [&lengths]( const api::PatchSide& side ) -> SmallVector<std::variant<bool, size_t>, 3> {
            switch( side )
            {
                case api::PatchSide::S0: return { size_t{ 0 }, true };
                case api::PatchSide::S1: return { lengths.at( 0 ) - 1, true };
                case api::PatchSide::T0: return { true, size_t{ 0 } };
                case api::PatchSide::T1: return { true, lengths.at( 1 ) - 1 };
                default: throw std::runtime_error( "Invalid side for 2D patch." );
            }
        }( side );

        std::vector<size_t> out;
        out.reserve( lengths.at( is_S ? 0 : 1 ) );

        const size_t offset = is_S ? component_bases.at( 0 )->numFunctions() : 0;

        util::iterateTensorProduct( lengths, { 0, 1 }, iter_dir, [&]( const util::IndexVec& iv ) {
            out.emplace_back( util::flatten( iv, lengths ) + offset );
        } );

        return out;
    };

    const std::vector<size_t> dirichlet_funcs = util::concatenate( {
        tangent_on_side( api::PatchSide::S0 ),
        tangent_on_side( api::PatchSide::S1 ),
        tangent_on_side( api::PatchSide::T0 ),
        tangent_on_side( api::PatchSide::T1 )
    } );

    Eigen::VectorXi id_array = Eigen::VectorXi::Zero( drcd.HCURL_ss.numFunctions() + drcd.H1_ss.numFunctions() );
    int j = 0;
    for( size_t i = 0; i < drcd.HCURL_ss.numFunctions() + drcd.H1_ss.numFunctions(); i++ )
    {
        if( std::find( dirichlet_funcs.begin(), dirichlet_funcs.end(), i ) != dirichlet_funcs.end() )//or i == drcd.HCURL_ss.numFunctions() )
        {
            id_array( i ) = -1;
        }
        else
        {
            id_array( i ) = j++;
        }
    }

    return id_array;
}

Eigen::VectorXi idArray( const DeRhamComplexHierarchicalDiscretization& drcd )
{
    const auto& hcurl = drcd.HCURL_ss;
    const auto& component_bases = hcurl.scalarBases();
    const size_t num_levels = component_bases.at( 0 )->basisComplex().parametricAtlas().cmap().numLevels();
    SmallVector<size_t, 4> offsets;
    offsets.push_back( 0 );
    for( size_t i = 0; i < component_bases.size(); i++ )
    {
        offsets.push_back( offsets.back() + component_bases.at( i )->numFunctions() );
    }
    std::cout << "offsets: " << offsets << std::endl;
    const util::IndexVec order = hcurl.basisComplex().parametricAtlas().cmap().dim() == 2 ? util::IndexVec{ 0, 1 } : util::IndexVec{ 0, 1, 2 };
    const auto get_iter_dir = []( const util::IndexVec& lengths, const api::PatchSide& side ) -> SmallVector<std::variant<bool, size_t>, 3> {
        if( lengths.size() == 2 )
        {
            switch( side )
            {
                case api::PatchSide::S0: return { size_t{ 0 }, true };
                case api::PatchSide::S1: return { lengths.at( 0 ) - 1, true };
                case api::PatchSide::T0: return { true, size_t{ 0 } };
                case api::PatchSide::T1: return { true, lengths.at( 1 ) - 1 };
                default: throw std::runtime_error( "Invalid side for 2D patch." );
            }
        }
        else if( lengths.size() == 3 )
        {
            switch( side )
            {
                case api::PatchSide::S0: return { size_t{ 0 }, true, true };
                case api::PatchSide::S1: return { lengths.at( 0 ) - 1, true, true };
                case api::PatchSide::T0: return { true, size_t{ 0 }, true };
                case api::PatchSide::T1: return { true, lengths.at( 1 ) - 1, true };
                case api::PatchSide::U0: return { true, true, size_t{ 0 } };
                case api::PatchSide::U1: return { true, true, lengths.at( 2 ) - 1 };
            }
        }
        else
            throw std::runtime_error( "Invalid number of dimensions for patch." );
    };
    const auto tangent_on_side = [&]( const api::PatchSide& side ) {
        const size_t which_coord = side == api::PatchSide::S0 or side == api::PatchSide::S1 ? 0 : side == api::PatchSide::T0 or side == api::PatchSide::T1 ? 1 : 2;

        std::vector<size_t> result;
        for( size_t comp_ii = 0; comp_ii < component_bases.size(); comp_ii++ )
        {
            if( comp_ii == which_coord )
                continue;

            const basis::HierarchicalTPSplineSpace& hier_comp = *component_bases.at( comp_ii );
            for( size_t i = 0; i < num_levels; i++ )
            {
                const basis::TPSplineSpace& comp = *hier_comp.refinementLevels().at( i );
                const util::IndexVec lengths = getTPLengths( comp );
                const auto iter_dir = get_iter_dir( lengths, side );
                const auto& exop = hier_comp.levelExtractionOperators().at( i );

                util::iterateTensorProduct( lengths, order, iter_dir, [&]( const util::IndexVec& iv ) {
                    const size_t fid = util::flatten( iv, lengths );
                    for( Eigen::SparseMatrix<double>::InnerIterator it( exop, fid ); it; ++it )
                    {
                        result.push_back( it.row() + offsets.at( comp_ii ) );
                        if( result.back() == 544 ) 
                        {
                            std::cout << "side: " << std::to_underlying( side ) << " comp_ii: " << comp_ii << " " << iter_dir << " " << iv << std::endl;
                            std::cout << "level: " << i << " lengths: " << lengths << std::endl;
                            std::cout << "fid: " << fid << std::endl;
                        }
                    }
                } );
            }
        }

        std::ranges::sort( result );
        result.erase( std::ranges::unique( result ).end(), result.end() );
        return result;
    };

    std::vector<size_t> h1_bcs;
    const auto h1_on_side = [&]( const api::PatchSide& side ) {
        std::vector<size_t> result;
        for( size_t i = 0; i < num_levels; i++ )
        {
            const basis::TPSplineSpace& comp = *drcd.H1_ss.refinementLevels().at( i );
            const util::IndexVec lengths = getTPLengths( comp );
            const auto iter_dir = get_iter_dir( lengths, side );
            const auto& exop = drcd.H1_ss.levelExtractionOperators().at( i );

            util::iterateTensorProduct( lengths, order, iter_dir, [&]( const util::IndexVec& iv ) {
                const size_t fid = util::flatten( iv, lengths );
                for( Eigen::SparseMatrix<double>::InnerIterator it( exop, fid ); it; ++it )
                {
                    result.push_back( it.row() + offsets.back() );
                }
            } );
        }
        std::ranges::sort( result );
        result.erase( std::ranges::unique( result ).end(), result.end() );
        return result;
    };

    const std::vector<size_t> dirichlet_funcs = component_bases.size() == 2 ? util::concatenate( {
        tangent_on_side( api::PatchSide::S0 ),
        tangent_on_side( api::PatchSide::S1 ),
        tangent_on_side( api::PatchSide::T0 ),
        tangent_on_side( api::PatchSide::T1 ),
        h1_on_side( api::PatchSide::S0 ),
        h1_on_side( api::PatchSide::S1 ),
        h1_on_side( api::PatchSide::T0 ),
        h1_on_side( api::PatchSide::T1 )
    } ) : util::concatenate( {
        tangent_on_side( api::PatchSide::S0 ),
        tangent_on_side( api::PatchSide::S1 ),
        tangent_on_side( api::PatchSide::T0 ),
        tangent_on_side( api::PatchSide::T1 ),
        tangent_on_side( api::PatchSide::U0 ),
        tangent_on_side( api::PatchSide::U1 ),
        h1_on_side( api::PatchSide::S0 ),
        h1_on_side( api::PatchSide::S1 ),
        h1_on_side( api::PatchSide::T0 ),
        h1_on_side( api::PatchSide::T1 ),
        h1_on_side( api::PatchSide::U0 ),
        h1_on_side( api::PatchSide::U1 )
    } );

    Eigen::VectorXi id_array = Eigen::VectorXi::Zero( drcd.HCURL_ss.numFunctions() + drcd.H1_ss.numFunctions() );
    int j = 0;
    for( size_t i = 0; i < drcd.HCURL_ss.numFunctions() + drcd.H1_ss.numFunctions(); i++ )
    {
        if( std::find( dirichlet_funcs.begin(), dirichlet_funcs.end(), i ) != dirichlet_funcs.end() )//or i == drcd.HCURL_ss.numFunctions() )
        {
            id_array( i ) = -1;
        }
        else
        {
            id_array( i ) = j++;
        }
    }

    return id_array;
}

int main(int argc, char** argv)
{
    PetscCall( SlepcInitialize(&argc, &argv, NULL, NULL) );
    // See arXiv:2411.15828v1 for the strong form of the Maxwell eigenproblem

    // DeRhamComplexHierarchicalDiscretization drcd = cases( TestCase::Case2p1KittyCornersIntersection_2d, 2 );
    DeRhamComplexHierarchicalDiscretization drcd = cases( TestCase::Case2d_nForm );
    // DeRhamComplexHierarchicalDiscretization drcd = cases( TestCase::ShortestChainExists_2d );

    const size_t dim = drcd.H1_ss.basisComplex().parametricAtlas().cmap().dim();
    const param::ParentDomain domain = param::cubeDomain( dim );
    const size_t degree = [&](){
        std::optional<topology::Cell> c;
        iterateCellsWhile( drcd.H1_ss.basisComplex().parametricAtlas().cmap(), drcd.H1_ss.basisComplex().parametricAtlas().cmap().dim(), [&]( const topology::Cell& cell ) {
            c = cell;
            return false;
        } );
        return drcd.H1_ss.basisComplex().parentBasis( c.value() ).mBasisGroups.at( 0 ).degrees.at( 0 );
    }();
    const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( degree + 1, domain );

    const MatrixX3dMax cpts_t = drcd.controlPoints().transpose();
    const io::BezierOutputObject bz_out( drcd.getH1().splineSpace(), cpts_t );
    io::outputBezierMeshToVTK( bz_out, "maxwell_eigenproblem_initial_mesh.vtu" );

    io::outputPartialBezierMeshToVTK( bz_out, "maxwell_eigenproblem_initial_mesh_partial.vtu", [&]( const std::function<void(const topology::Cell& )>& callback ) {
        iterateCellsWhile( drcd.H1_ss.basisComplex().parametricAtlas().cmap(), dim, [&]( const topology::Cell& cell ) {
            const size_t level = drcd.H1_ss.basisComplex().parametricAtlas().cmap().unrefinedAncestorDartOfCell( cell ).first;
            if( level == 1 )
                callback( cell );
            return true;
        } );
    } );

    const Eigen::VectorXi id_array = idArray( drcd );
    std::cout << "ID array:\n";// << id_array.transpose() << std::endl << std::endl;
    // const Eigen::VectorXi id_array = Eigen::VectorXi::LinSpaced( drcd.getHCURL().splineSpace().numFunctions() + drcd.getH1().splineSpace().numFunctions(), 0, drcd.getHCURL().splineSpace().numFunctions() + drcd.getH1().splineSpace().numFunctions() - 1 );

    const size_t n_hcurl_constraints = ( id_array.head( drcd.HCURL_ss.numFunctions() ).array() == -1 ).count();
    const size_t n_constraints = ( id_array.array() == -1 ).count();
    const size_t n_h1_constraints = n_constraints - n_hcurl_constraints;
    std::cout << drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints << " free HCURL functions\n";
    std::cout << drcd.getH1().splineSpace().numFunctions() - n_h1_constraints << " free H1 functions\n";

    const Eigen::SparseMatrix<double> K = assembleStiffnessMatrix( drcd, quad_rule, id_array );
    const Eigen::SparseMatrix<double> M = massMatrix( drcd, quad_rule, id_array );
    std::cout << "Assembled stiffness and mass matrices." << std::endl;

    Eigen::saveMarket( K, "stiffness_matrix.mtx" );
    Eigen::saveMarket( M, "mass_matrix.mtx" );

    // std::vector<basis::FunctionId> bc_funcs;
    // for( Eigen::Index i = 0; i < id_array.size(); i++ )
    // {
    //     if( id_array( i ) == -1 )
    //     {
    //         bc_funcs.push_back( basis::FunctionId( i ) );
    //     }
    // }
    // std::cout << bc_funcs << std::endl;
    // io::outputVectorBasis( drcd.getHCURL().splineSpace(), drcd.H1_ss, drcd.controlPoints(), "test", bc_funcs );

    // std::cout << K.toDense() << std::endl;

    // std::cout << "Stiffness matrix K:" << std::endl;
    // std::cout << Eigen::MatrixXd( K ) << std::endl;
    // std::cout << "Mass matrix M:" << std::endl;
    // std::cout << Eigen::MatrixXd( M ) << std::endl;
    // std::cout << "=============================================" << std::endl;
    
    // // // Solve generalized eigenproblem
    // Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver( K.toDense(), M.toDense() );

    // // std::cout << "Eigenvalues:" << std::endl;
    // // std::cout << solver.eigenvalues() << std::endl;
    // std::cout << solver.eigenvalues().transpose() << std::endl;
    // Eigen::VectorXd magnitudes = solver.eigenvalues().cwiseAbs();
    // std::sort( magnitudes.data(), magnitudes.data() + magnitudes.size() );
    // std::cout << "Eigenvalue magnitudes:" << std::endl;
    // std::cout << magnitudes.tail( magnitudes.size() - (magnitudes.array() <= 1e-10 ).count()).transpose() << std::endl;

    // std::cout << (magnitudes.array() <= 1e-10 ).count() << " zero eigenvalues due to constraints." << std::endl;
    // std::cout << solver.alphas().transpose() << std::endl;
    // std::cout << solver.betas().transpose() << std::endl;
    // std::cout << "Total eigenvalues: " << magnitudes.size() << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "=============================================" << std::endl;
    // std::cout << "Constraints: " << n_constraints << std::endl;
    // std::cout << id_array.transpose() << std::endl;
    std::cout << "=============================================" << std::endl;
    // std::cout << "HCURL block size: " << drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints << std::endl;

    // const Eigen::MatrixXd k_part = K.toDense().topLeftCorner( drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints,
    //                                                                     drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints );
    // const Eigen::MatrixXd m_part = M.toDense().topLeftCorner( drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints,
    //                                                                     drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints );

    // // std::cout << "K (HCURL block only):" << std::endl;
    // // std::cout << k_part << std::endl;
    // // std::cout << "M (HCURL block only):" << std::endl;
    // // std::cout << m_part << std::endl;

    // Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver2( k_part, m_part );
    // // std::cout << "Eigenvalues (HCURL block only):" << std::endl;
    // // std::cout << solver2.eigenvalues().transpose() << std::endl;
    // Eigen::VectorXd magnitudes2 = solver2.eigenvalues().cwiseAbs();
    // std::sort( magnitudes2.data(), magnitudes2.data() + magnitudes2.size() );
    // std::cout << "Eigenvalue magnitudes (HCURL block only):" << std::endl;
    // std::cout << magnitudes2.tail( magnitudes2.size() - (magnitudes2.array() <= 1e-10 ).count()).transpose() << std::endl;

    // std::cout << (magnitudes2.array() <= 1e-10 ).count() << " zero eigenvalues due to constraints." << std::endl;
    // std::cout << "Total eigenvalues: " << magnitudes2.size() << std::endl;

    // Eigen::Index idx;
    // (solver2.eigenvalues().cwiseAbs().array() - magnitudes2((magnitudes2.array() <= 1e-10 ).count()+3)).abs().matrix().minCoeff( &idx );
    // const Eigen::VectorXd omega_known = solver2.eigenvectors().col( idx ).real();

    // // const Eigen::VectorXd omega_analytical = calculateDivFreeField( drcd, id_array );

    // const auto constraint_residual = K.toDense().bottomLeftCorner( drcd.H1_ss.numFunctions() - 1, drcd.HCURL_ss.numFunctions() - n_hcurl_constraints ) * omega_known;
    // std::cout << "Constraint residual (should be close to 0):\n" << constraint_residual.norm() << std::endl;


    try {
        // Solve eigenvalue problem
        const size_t nev = 200;//drcd.getHCURL().splineSpace().numFunctions() / 2;
        std::cout << "Solving generalized eigenvalue problem for " << nev << " eigenvalues" << std::endl;

        auto [eigenvalues, eigenvectors] = slepc_utils::solveGeneralizedEigenvalueSparse(K, M, nev, drcd.HCURL_ss.numFunctions() - n_hcurl_constraints );

        // Print results
        std::cout << "\nEigenvalues:" << std::endl;
        std::sort( eigenvalues.data(), eigenvalues.data() + eigenvalues.size() );
        const Eigen::Map<Eigen::VectorXd> eigenvalues_sorted( eigenvalues.data(), eigenvalues.size() );
        const size_t num_zero_eigenvalues = (eigenvalues_sorted.array() <= 1e-10 ).count();
        std::cout << num_zero_eigenvalues << " zero eigenvalues due to constraints." << std::endl;
        std::cout << "Total eigenvalues: " << eigenvalues.size() << std::endl;
        std::cout << "Non-zero eigenvalue magnitudes:" << std::endl;
        std::cout << Eigen::Map<Eigen::VectorXd>( eigenvalues.data(), eigenvalues.size() ).tail( eigenvalues.size() - num_zero_eigenvalues ).transpose() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        PetscCall( SlepcFinalize() );
        return 1;
    }

    // const Eigen::SparseMatrix<double> k_part = K.topLeftCorner( drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints,
    //                                                                     drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints );
    // const Eigen::SparseMatrix<double> m_part = M.topLeftCorner( drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints,
    //                                                                     drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints );
    // const Eigen::SparseMatrix<double> b_part = K.bottomLeftCorner( drcd.getH1().splineSpace().numFunctions(),
    //                                                                     drcd.getHCURL().splineSpace().numFunctions() - n_hcurl_constraints );

    // return assembly::solveMixedEigenproblem( argc, argv, 100, k_part, m_part, b_part );
    PetscCall( SlepcFinalize() );

    return 0;
}
