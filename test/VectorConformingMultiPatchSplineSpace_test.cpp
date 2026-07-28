#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <VectorConformingMultiPatchSplineSpace.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CommonUtils.hpp>
#include <numbers>
#include <Eigen/Dense>

using namespace topology;
using namespace param;
using namespace basis;

Eigen::MatrixXd getTransformedEvals( const eval::SplineSpaceEvaluator& evaler,
                                     const eval::SplineSpaceEvaluator& geom_evaler,
                                     const Eigen::MatrixXd& geom,
                                     const ConformingType conforming_type )
{
    if( conforming_type == ConformingType::Divergence )
    {
        return eval::evaluateSpatialHDivBasisValues( evaler, geom_evaler, geom );
    }
    else if( conforming_type == ConformingType::Curl )
    {
        return eval::evaluateSpatialHCurlBasisValues( evaler, geom_evaler, geom );
    }
    else
    {
        throw std::runtime_error( "getTransformedEvals: Unsupported conforming type!" );
    }
}

void test_c0( const VectorConformingMultiPatchSplineSpace& ss,
              const MultiPatchSplineSpace& geom_ss,
              const Eigen::MatrixXd& geom,
              const size_t n_funcs_expected,
              const size_t num_overlapping_expected,
              const ConformingType conforming_type )
{
    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto& cmap = static_cast<const MultiPatchCombinatorialMap&>( ss.basisComplex().parametricAtlas().cmap() );
    CHECK( ss.numFunctions() == n_funcs_expected );
    CHECK( ss.numVectorComponents() == cmap.dim() );
    
    // Check functions are connected
    eval::SplineSpaceEvaluator evaler( ss, 0 );
    eval::SplineSpaceEvaluator geom_evaler( geom_ss, 1 );
    for( size_t i = 0; i < cmap.constituents().size(); i++ )
    {
        const auto& submap = cmap.constituents().at( i );
        iterateCellsWhile( *submap, cmap.dim() - 1, [&]( const Cell& e ) {
            if( onBoundary( *submap, e.dart() ) )
            {
                const Dart one_side_d = cmap.toGlobalDart( i, e.dart() );
                const auto maybe_phi2 = phi( cmap, cmap.dim(), one_side_d );
                if( maybe_phi2.has_value() )
                {
                    // Evaluate on both sides of the edge
                    const Cell one_side_f( one_side_d, cmap.dim() );
                    const Cell one_side_e( one_side_d, cmap.dim() - 1 );
                    const ParentPoint one_side_ppt = pointOnBoundary( cubeDomain( cmap.dim() ), parentDomainBoundary( atlas, one_side_e ) );

                    evaler.localizeElement( one_side_f );
                    evaler.localizeParentPoint( one_side_ppt );
                    geom_evaler.localizeElement( one_side_f );
                    geom_evaler.localizeParentPoint( one_side_ppt );
                    const Eigen::MatrixXd one_side_evals = getTransformedEvals( evaler, geom_evaler, geom, conforming_type );
                    const auto one_side_conn = ss.connectivity( one_side_f );

                    const Cell other_side_f( maybe_phi2.value(), cmap.dim() );
                    const Cell other_side_e( maybe_phi2.value(), cmap.dim() - 1 );
                    const ParentPoint other_side_ppt = pointOnBoundary( cubeDomain( cmap.dim() ), parentDomainBoundary( atlas, other_side_e ) );

                    evaler.localizeElement( other_side_f );
                    evaler.localizeParentPoint( other_side_ppt );
                    geom_evaler.localizeElement( other_side_f );
                    geom_evaler.localizeParentPoint( other_side_ppt );
                    const Eigen::MatrixXd other_side_evals = getTransformedEvals( evaler, geom_evaler, geom, conforming_type );
                    const auto other_side_conn = ss.connectivity( other_side_f );

                    size_t n_overlapping_funcs = 0;
                    for( size_t i = 0; i < one_side_conn.size(); i++ )
                    {
                        for( size_t j = 0; j < other_side_conn.size(); j++ )
                        {
                            if( one_side_conn.at( i ) == other_side_conn.at( j ) )
                            {
                                CHECK( util::equals( one_side_evals.row( i ), other_side_evals.row( j ), 1e-9 ) );
                                n_overlapping_funcs++;
                            }
                        }
                    }
                    CHECK( n_overlapping_funcs == num_overlapping_expected );
                }
            }
            return true;
        } );
    }
}

TEST_CASE( "Simple 2d two-patch divergence-conforming spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv2( { 0, 0, 0, 0, 1, 2, 2, 2, 2 }, ptol );
    const size_t degree1 = 2;
    const size_t degree2 = 3;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2}, {degree1, degree2} ) );
    const auto ss_tp_2 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv1}, {degree2, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& bc_tp_2 = ss_tp_2->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& atlas_tp_2 = bc_tp_2->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();
    const auto& cmap_tp_2 = atlas_tp_2->cmapPtr();

    const auto dc_bc_1 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_1, ConformingType::Divergence );
    const auto dc_bc_2 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_2, ConformingType::Divergence );

    const auto dc_ss_1 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_1, *ss_tp_1 );
    const auto dc_ss_2 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_2, *ss_tp_2 );

    SECTION( "Positive S with Positive T" )
    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 4 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
            const Eigen::MatrixXd patch2_geom =
                ( util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                  Eigen::Rotation2Dd( -1.0 * std::numbers::pi / 2.0 ).matrix().transpose() )
                    .rowwise() +
                Eigen::RowVector2d( 1.0, 2.0 );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Divergence );

        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );
        test_c0( ss, geom_ss, geom, 40, degree2, ConformingType::Divergence );
    }

    SECTION( "Positive S with Negative T" )
    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 2 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
            const Eigen::MatrixXd patch2_geom =
                ( util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                  Eigen::Rotation2Dd( std::numbers::pi / 2.0 ).matrix().transpose() )
                    .rowwise() +
                Eigen::RowVector2d( 0.0, 2.0 );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Divergence );

        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );
        test_c0( ss, geom_ss, geom, 40, degree2, ConformingType::Divergence );
    }
}

TEST_CASE( "Simple 2d two-patch curl-conforming spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv2( { 0, 0, 0, 0, 1, 2, 2, 2, 2 }, ptol );
    const size_t degree1 = 2;
    const size_t degree2 = 3;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2}, {degree1, degree2} ) );
    const auto ss_tp_2 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv1}, {degree2, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& bc_tp_2 = ss_tp_2->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& atlas_tp_2 = bc_tp_2->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();
    const auto& cmap_tp_2 = atlas_tp_2->cmapPtr();

    const auto dc_bc_1 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_1, ConformingType::Curl );
    const auto dc_bc_2 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_2, ConformingType::Curl );

    const auto dc_ss_1 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_1, *ss_tp_1 );
    const auto dc_ss_2 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_2, *ss_tp_2 );

    SECTION( "Positive T with Negative S" )
    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 4 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
            const Eigen::MatrixXd patch2_geom =
                ( util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                  Eigen::Rotation2Dd( -1.0 * std::numbers::pi / 2.0 ).matrix().transpose() )
                    .rowwise() +
                Eigen::RowVector2d( 1.0, 2.0 );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Curl );

        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );
        test_c0( ss, geom_ss, geom, 40, degree2, ConformingType::Curl );
    }

    SECTION( "Positive T with Positive S" )
    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 2 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
            const Eigen::MatrixXd patch2_geom =
                ( util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                  Eigen::Rotation2Dd( std::numbers::pi / 2.0 ).matrix().transpose() )
                    .rowwise() +
                Eigen::RowVector2d( 2.0, 0.0 );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Curl );

        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );
        test_c0( ss, geom_ss, geom, 40, degree2, ConformingType::Curl );
    }
}

TEST_CASE( "Simple 3d two-patch div conforming spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv2( { 0, 0, 0, 0, 1, 2, 2, 2, 2 }, ptol );
    const size_t degree1 = 2;
    const size_t degree2 = 3;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv1, kv1}, {degree2, degree1, degree1} ) );
    const auto ss_tp_2 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2, kv1}, {degree1, degree2, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& bc_tp_2 = ss_tp_2->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& atlas_tp_2 = bc_tp_2->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();
    const auto& cmap_tp_2 = atlas_tp_2->cmapPtr();

    const auto dc_bc_1 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_1, ConformingType::Divergence );
    const auto dc_bc_2 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_2, ConformingType::Divergence );

    const auto dc_ss_1 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_1, *ss_tp_1 );
    const auto dc_ss_2 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_2, *ss_tp_2 );

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 19 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom =
                util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ) } );
            const Eigen::MatrixXd patch2_geom = util::tensorProduct(
                { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) * 
                ( Eigen::AngleAxisd( std::numbers::pi / 2.0, Eigen::Vector3d::UnitZ() ).toRotationMatrix() );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Divergence );
        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );

        const size_t num_functions = dc_ss_1->numFunctions() + dc_ss_2->numFunctions() - 8;
        test_c0( ss, geom_ss, geom, num_functions, degree2 * degree1, ConformingType::Divergence );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 6 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom =
                util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ) } );
            const Eigen::MatrixXd patch2_geom = (util::tensorProduct(
                { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                Eigen::AngleAxisd( std::numbers::pi / 2.0, Eigen::Vector3d::UnitZ() ).toRotationMatrix() *
                 Eigen::AngleAxisd( -1 * std::numbers::pi / 2.0, Eigen::Vector3d::UnitX() ).toRotationMatrix() ).rowwise() + Eigen::RowVector3d( 0.0, 0.0, 1.0 );

            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Divergence );
        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );

        const size_t num_functions = dc_ss_1->numFunctions() + dc_ss_2->numFunctions() - 8;
        test_c0( ss, geom_ss, geom, num_functions, degree2 * degree1, ConformingType::Divergence );
    }
}

TEST_CASE( "Simple 3d two-patch curl conforming spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 0, 1, 1, 1 }, ptol );
    const KnotVector kv2( { 0, 0, 0, 0, 1, 2, 2, 2, 2 }, ptol );
    const size_t degree1 = 2;
    const size_t degree2 = 3;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv1, kv1}, {degree2, degree1, degree1} ) );
    const auto ss_tp_2 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2, kv1}, {degree1, degree2, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& bc_tp_2 = ss_tp_2->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& atlas_tp_2 = bc_tp_2->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();
    const auto& cmap_tp_2 = atlas_tp_2->cmapPtr();

    const auto dc_bc_1 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_1, ConformingType::Curl );
    const auto dc_bc_2 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_2, ConformingType::Curl );

    const auto dc_ss_1 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_1, *ss_tp_1 );
    const auto dc_ss_2 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_2, *ss_tp_2 );

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 19 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom =
                util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ) } );
            const Eigen::MatrixXd patch2_geom = util::tensorProduct(
                { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) * 
                ( Eigen::AngleAxisd( std::numbers::pi / 2.0, Eigen::Vector3d::UnitZ() ).toRotationMatrix() );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Curl );
        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );

        const size_t num_functions = dc_ss_1->numFunctions() + dc_ss_2->numFunctions() - 22;
        test_c0( ss, geom_ss, geom, num_functions, 17, ConformingType::Curl );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 6 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_2 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom =
                util::tensorProduct( { grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ) } );
            const Eigen::MatrixXd patch2_geom = (util::tensorProduct(
                { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv1, degree1 ) } ) *
                Eigen::AngleAxisd( std::numbers::pi / 2.0, Eigen::Vector3d::UnitZ() ).toRotationMatrix() *
                 Eigen::AngleAxisd( -1 * std::numbers::pi / 2.0, Eigen::Vector3d::UnitX() ).toRotationMatrix() ).rowwise() + Eigen::RowVector3d( 0.0, 0.0, 1.0 );

            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Curl );
        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_2 } );

        const size_t num_functions = dc_ss_1->numFunctions() + dc_ss_2->numFunctions() - 22;
        test_c0( ss, geom_ss, geom, num_functions, 17, ConformingType::Curl );
    }
}

TEST_CASE( "3d three-patch curl conforming spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 0, 1, 1, 1 }, ptol );
    const size_t degree1 = 2;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv1, kv1}, {degree1, degree1, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();

    const auto dc_bc_1 = std::make_shared<const VectorConformingBasisComplex>( bc_tp_1, ConformingType::Curl );

    const auto dc_ss_1 = std::make_shared<const VectorConformingTPSplineSpace>( dc_bc_1, *ss_tp_1 );

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_1, cmap_tp_1 }, { { { 0, Dart( 13 ) }, { 1, Dart( 1 ) } }, { { 0, Dart( 7 ) }, { 2, Dart( 19 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_1, atlas_tp_1 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_1, bc_tp_1 } );

        const MultiPatchSplineSpace geom_ss( bc, { ss_tp_1, ss_tp_1, ss_tp_1 } );

        const Eigen::MatrixXd geom = [&]() -> Eigen::MatrixXd {
            const Eigen::MatrixXd patch1_geom =
                util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ), grevillePoints( kv1, degree1 ) } );
            const Eigen::MatrixXd patch2_geom = patch1_geom.rowwise() + Eigen::RowVector3d( 0.0, 1.0, 0.0 );
            const Eigen::MatrixXd patch3_geom = patch1_geom.rowwise() + Eigen::RowVector3d( 1.0, 0.0, 0.0 );
            return multiPatchCoefficients( geom_ss, {patch1_geom, patch2_geom, patch3_geom} ).transpose();
        }();

        const auto dc_bc = std::make_shared<const VectorConformingBasisComplex>( bc, ConformingType::Curl );
        const VectorConformingMultiPatchSplineSpace ss( dc_bc, { dc_ss_1, dc_ss_1, dc_ss_1 } );

        const size_t num_functions = 3 * dc_ss_1->numFunctions() - 24;
        test_c0( ss, geom_ss, geom, num_functions, 12, ConformingType::Curl );
    }
}
