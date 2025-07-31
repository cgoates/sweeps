#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CommonUtils.hpp>
#include <MFEMOutput.hpp>
#include <fstream>

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

void test_c0( const MultiPatchSplineSpace& ss,
              const size_t n_funcs_expected,
              const size_t num_overlapping_expected )
{
    const auto& cmap = ss.basisComplex().parametricAtlas().cmap();

    for( const auto& submap : cmap.constituents() )
        test_c0( ss, *submap, n_funcs_expected, num_overlapping_expected );
}

TEST_CASE( "Simple 2d two-patch spline space" )
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

        const MultiPatchSplineSpace ss( bc, { ss_tp_1, ss_tp_2 } );

        test_c0( ss, *cmap_tp_1, 25, degree2 + 1 );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 2 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss( bc, { ss_tp_1, ss_tp_2 } );

        test_c0( ss, *cmap_tp_1, 25, degree2 + 1 );
    }
}

TEST_CASE( "Simple 3d two-patch spline space" )
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

        const MultiPatchSplineSpace ss( bc, { ss_tp_1, ss_tp_2 } );

        test_c0( ss, *cmap_tp_1, 75, ( degree2 + 1 ) * ( degree1 + 1 ) );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 1 ) }, { 1, Dart( 9 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

        const MultiPatchSplineSpace ss( bc, { ss_tp_1, ss_tp_2 } );

        test_c0( ss, *cmap_tp_1, 75, ( degree2 + 1 ) * ( degree1 + 1 ) );
    }
}

TEST_CASE( "MFEM Output" )
{
    const double ptol = 1e-9;
    const KnotVector kv1( { 0, 0, 1, 1 }, ptol );
    const auto& kv2 = kv1;
    const size_t degree1 = 1;
    const size_t degree2 = degree1;

    const auto ss_tp_1 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv1, kv1}, {degree2, degree1, degree1} ) );
    const auto ss_tp_2 = std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2, kv1}, {degree1, degree2, degree1} ) );
    const auto& bc_tp_1 = ss_tp_1->basisComplexPtr();
    const auto& bc_tp_2 = ss_tp_2->basisComplexPtr();
    const auto& atlas_tp_1 = bc_tp_1->parametricAtlasPtr();
    const auto& atlas_tp_2 = bc_tp_2->parametricAtlasPtr();
    const auto& cmap_tp_1 = atlas_tp_1->cmapPtr();
    const auto& cmap_tp_2 = atlas_tp_2->cmapPtr();

    const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
        MultiPatchCombinatorialMap( { cmap_tp_1, cmap_tp_2 }, { { { 0, Dart( 7 ) }, { 1, Dart( 19 ) } } } ) );
    const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
        cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp_1, atlas_tp_2 } );
    const auto bc = std::make_shared<const MultiPatchBasisComplex>(
        atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp_1, bc_tp_2 } );

    const MultiPatchSplineSpace ss( bc, { ss_tp_1, ss_tp_2 } );
    
    const Eigen::MatrixXd cpts0 = ( Eigen::MatrixXd( 8, 3 ) <<
        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        1, 1, 0,
        0, 0, 1,
        1, 0, 1,
        0, 1, 1,
        1, 1, 1 ).finished();
    const Eigen::MatrixXd cpts1 = ( Eigen::MatrixXd( 8, 3 ) <<
        1, 0, 0,
        2, 0, 0,
        1, 1, 0,
        2, 1, 0,
        1, 0, 1,
        2, 0, 1,
        1, 1, 1,
        2, 1, 1 ).finished();
    
    const Eigen::MatrixXd cpts = multiPatchCoefficients( ss, { cpts0, cpts1 } );

    const std::string filename = "simple_3d_two_patch.mesh";

    io::outputMultiPatchSplinesToMFEM( ss, cpts.transpose(), filename );

    std::ifstream file( filename );
    if( not file.is_open() )
    {
        FAIL( "Could not open file " + filename );
    }

    std::string line;
    int patches = 0;
    
    // Read dimension from line 4
    for( int i = 1; i <= 4; ++i )
    {
        if( !std::getline( file, line ) )
        {
            FAIL( "File too short, couldn't read line 4 to find dimension." );
        }
    }
    CHECK( std::stoi( line ) == 3 );

    while( std::getline( file, line ) )
    {
        line.erase( 0, line.find_first_not_of( " \t" ) );
        line.erase( line.find_last_not_of( " \t" ) + 1 );

        if( line == "elements" )
        {
            std::getline( file, line );
            CHECK( std::stoi( line ) == 2 );
        }
        else if( line == "boundary" )
        {
            std::getline( file, line );
            CHECK( std::stoi( line ) == 10 );
        }
        else if( line == "edges" )
        {
            std::getline( file, line );
            CHECK( std::stoi( line ) == 20 );
        }
        else if( line == "vertices" )
        {
            std::getline( file, line );
            CHECK( std::stoi( line ) == 12 );
        }
        else if( line == "knotvectors" )
        {
            patches++;
        }
    }

    CHECK( patches == 2 );

    file.close();
}

TEST_CASE( "More complex 2d multipatch spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv( { 0, 0, 0, 1, 2, 3, 3, 3 }, ptol );
    const size_t degree = 2;

    const auto ss_tp = std::make_shared<const TPSplineSpace>( buildBSpline( {kv, kv}, {degree, degree} ) );
    const auto& bc_tp = ss_tp->basisComplexPtr();
    const auto& atlas_tp = bc_tp->parametricAtlasPtr();
    const auto& cmap_tp = atlas_tp->cmapPtr();

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp, cmap_tp, cmap_tp },
                                        { { { 0, Dart( 26 ) }, { 1, Dart( 0 ) } },
                                          { { 0, Dart( 27 ) }, { 2, Dart( 8 ) } },
                                          { { 1, Dart( 3 ) }, { 2, Dart( 9 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp, atlas_tp, atlas_tp } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp, bc_tp, bc_tp } );

        const MultiPatchSplineSpace ss( bc, { ss_tp, ss_tp, ss_tp } );

        test_c0( ss, 61, degree + 1 );
    }

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( std::vector<std::shared_ptr<const topology::TPCombinatorialMap>>( 8, cmap_tp ),
                                        { { { 0, Dart( 9 ) }, { 1, Dart( 3 ) } },
                                          { { 0, Dart( 26 ) }, { 3, Dart( 0 ) } },
                                          { { 1, Dart( 9 ) }, { 2, Dart( 3 ) } },
                                          { { 2, Dart( 26 ) }, { 4, Dart( 0 ) } },
                                          { { 3, Dart( 26 ) }, { 5, Dart( 0 ) } },
                                          { { 4, Dart( 26 ) }, { 7, Dart( 0 ) } },
                                          { { 5, Dart( 9 ) }, { 6, Dart( 3 ) } },
                                          { { 6, Dart( 9 ) }, { 7, Dart( 3 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>( 8, atlas_tp ) );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>( 8, bc_tp ) );

        const MultiPatchSplineSpace ss( bc, std::vector<std::shared_ptr<const TPSplineSpace>>( 8, ss_tp ) );

        test_c0( ss, 160, degree + 1 );
    }
}

TEST_CASE( "More complex 3d multipatch spline space" )
{
    const double ptol = 1e-9;
    const KnotVector kv( { 0, 0, 0, 1, 1, 1 }, ptol );
    const size_t degree = 2;

    const auto ss_tp = std::make_shared<const TPSplineSpace>( buildBSpline( {kv, kv, kv}, {degree, degree, degree} ) );
    const auto& bc_tp = ss_tp->basisComplexPtr();
    const auto& atlas_tp = bc_tp->parametricAtlasPtr();
    const auto& cmap_tp = atlas_tp->cmapPtr();

    {
        const auto cmap = std::make_shared<const MultiPatchCombinatorialMap>(
            MultiPatchCombinatorialMap( { cmap_tp, cmap_tp, cmap_tp, cmap_tp },
                                        { { { 0, Dart( 0 ) }, { 1, Dart( 5 ) } },
                                          { { 0, Dart( 8 ) }, { 2, Dart( 5 ) } },
                                          { { 0, Dart( 16 ) }, { 3, Dart( 23 ) } },
                                          { { 1, Dart( 8 ) }, { 2, Dart( 22 ) } },
                                          { { 1, Dart( 16 ) }, { 3, Dart( 2 ) } },
                                          { { 2, Dart( 16 ) }, { 3, Dart( 8 ) } } } ) );
        const auto atlas = std::make_shared<const MultiPatchParametricAtlas>(
            cmap, std::vector<std::shared_ptr<const TPParametricAtlas>>{ atlas_tp, atlas_tp, atlas_tp, atlas_tp } );
        const auto bc = std::make_shared<const MultiPatchBasisComplex>(
            atlas, std::vector<std::shared_ptr<const TPBasisComplex>>{ bc_tp, bc_tp, bc_tp, bc_tp } );

        const MultiPatchSplineSpace ss( bc, { ss_tp, ss_tp, ss_tp, ss_tp } );

        test_c0( ss, 65, ( degree + 1 ) * ( degree + 1 ) );
    }
}

TEST_CASE( "Single patch periodic with multipatch" )
{
    const size_t n_elems_s = 2;
    const size_t p = 2;
    const KnotVector kv = integerKnotsWithNElems( 1, p );
    const KnotVector kv2 = integerKnotsWithNElems( n_elems_s, p );
    const auto ss_3d = std::make_shared<const TPSplineSpace>( buildBSpline( {kv2, kv, kv}, {p, p, p} ) );

    const MultiPatchSplineSpace ss = buildMultiPatchSplineSpace( { ss_3d }, { { { 0, Dart( 3 * 6 + 1 ) }, { 0, Dart( ( 4 * ( n_elems_s - 1 ) + 1 ) * 6 + 1 ) } } } );
    CHECK( ss.numFunctions() == 27 );
}
