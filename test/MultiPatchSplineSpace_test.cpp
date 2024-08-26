#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CommonUtils.hpp>

using namespace topology;
using namespace param;
using namespace basis;

void test_c0( const MultiPatchSplineSpace& ss,
              const CombinatorialMap& submap,
              const size_t n_funcs_expected,
              const size_t num_overlapping_expected )
{
    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto& cmap = atlas.cmap();
    CHECK( ss.numFunctions() == n_funcs_expected );
    CHECK( ss.numVectorComponents() == 1 );
    
    // Check functions are connected
    eval::SplineSpaceEvaluator evaler( ss, 0 );
    iterateCellsWhile( submap, cmap.dim() - 1, [&]( const Cell& e ) {
        if( onBoundary( submap, e.dart() ) )
        {
            const Dart one_side_d = cmap.toGlobalDart( 0, e.dart() );
            const auto maybe_phi2 = phi( cmap, cmap.dim(), one_side_d );
            if( maybe_phi2.has_value() )
            {
                // Evaluate on both sides of the edge
                const Cell one_side_f( one_side_d, cmap.dim() );
                const Cell one_side_e( one_side_d, cmap.dim() - 1 );
                const ParentPoint one_side_ppt = pointOnBoundary( cubeDomain( cmap.dim() ), parentDomainBoundary( atlas, one_side_e ) );

                evaler.localizeElement( one_side_f );
                evaler.localizePoint( one_side_ppt );
                const Eigen::MatrixXd one_side_evals = evaler.evaluateBasis();
                const auto one_side_conn = ss.connectivity( one_side_f );

                const Cell other_side_f( maybe_phi2.value(), cmap.dim() );
                const Cell other_side_e( maybe_phi2.value(), cmap.dim() - 1 );
                const ParentPoint other_side_ppt = pointOnBoundary( cubeDomain( cmap.dim() ), parentDomainBoundary( atlas, other_side_e ) );

                evaler.localizeElement( other_side_f );
                evaler.localizePoint( other_side_ppt );
                const Eigen::MatrixXd other_side_evals = evaler.evaluateBasis();
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

TEST_CASE( "Simple 2d multipatch spline space" )
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

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 4 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss( bc,
                                        { ss_tp_1, ss_tp_2 },
                                        { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 },
                                          { 14, 11, 8, 5, 2, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 } } );

        test_c0( ss, *cmap_tp_1, 25, degree2 + 1 );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 2 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss( bc,
                                        { ss_tp_1, ss_tp_2 },
                                        { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 },
                                          { 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 2, 5, 8, 11, 14 } } );

        test_c0( ss, *cmap_tp_1, 25, degree2 + 1 );
    }
}

TEST_CASE( "Simple 3d multipatch spline space" )
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

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 19 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss(
            bc,
            { ss_tp_1, ss_tp_2 },
            { { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44 },
              { 0,  45, 46, 1,  47, 48, 2,  49, 50, 3,  51, 52, 4,  53, 54, 15, 55, 56, 16, 57, 58, 17, 59,
                60, 18, 61, 62, 19, 63, 64, 30, 65, 66, 31, 67, 68, 32, 69, 70, 33, 71, 72, 34, 73, 74 } } );

        test_c0( ss, *cmap_tp_1, 75, ( degree2 + 1 ) * ( degree1 + 1 ) );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 9 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss(
            bc,
            { ss_tp_1, ss_tp_2 },
            { { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44 },
              { 45,  46,  30,  47,  48,  31,  49,  50,  32,  51,  52, 33, 53, 54, 34, 55, 56, 15, 57, 58, 16, 59, 60,
                17, 61, 62, 18, 63, 64, 19, 65, 66, 0, 67, 68, 1, 69, 70, 2, 71, 72, 3, 73, 74, 4 } } );

        test_c0( ss, *cmap_tp_1, 75, ( degree2 + 1 ) * ( degree1 + 1 ) );
    }
}