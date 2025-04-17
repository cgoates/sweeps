#include <catch2/catch_test_macros.hpp>
#include <HierarchicalMultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <VTKOutput.hpp>
#include <SplineSpaceEvaluator.hpp>

using namespace topology;
using namespace param;
using namespace basis;

constexpr bool VTK_OUTPUT = false;

TEST_CASE( "Simple 2d multi patch hierarchical spline space 1" )
{
    const KnotVector kv1( {0, 0, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0.5, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0.25, 0.5, 0.75, 1, 1}, 1e-10 );

    const size_t degree = 1;
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv1, kv2 }, { degree, degree } ) );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv2, kv4 }, { degree, degree } ) );

    const auto mp_ss_1 = std::make_shared<const MultiPatchSplineSpace>( buildMultiPatchSplineSpace( { tp_ss_1, tp_ss_1 }, {
        { { 0, Dart( 1 ) }, { 1, Dart( 3 ) } }
    } ) );

    const auto mp_ss_2 = std::make_shared<const MultiPatchSplineSpace>( buildMultiPatchSplineSpace( { tp_ss_2, tp_ss_2 }, {
        { { 0, Dart( 1 ) }, { 1, Dart( 3 ) } }
    } ) );

    const HierarchicalMultiPatchSplineSpace ss = buildHierarchicalSplineSpace( { mp_ss_1, mp_ss_2 }, {
        { Face( 0 ), Face( 4 ), Face( 8 ) },
        { Face( 48 ), Face( 52 ), Face( 56 ), Face( 60 ) }
    } );

    CHECK( cellCount( ss.basisComplex().parametricAtlas().cmap(), 2 ) == 7 );
    CHECK( ss.numFunctions() == 12 );
    iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const Face& f ) {
        CHECK( ss.connectivity( f ).size() == 4 );
        return true;
    } );

    const Eigen::MatrixX2d geom = [&]() {
        const Eigen::MatrixXd greville_1 = grevillePoints( *ss.constituents().at( 0 ) );
        const Eigen::MatrixXd greville_2 = grevillePoints( *ss.constituents().at( 1 ) );

        Eigen::MatrixXd geom( ss.numFunctions(), 2 );

        for( size_t i = 0; i < ss.functionIdMap().at( 0 ).size(); i++ )
        {
            geom.row( ss.functionIdMap().at( 0 ).at( i ) ) = greville_1.row( i );
        }
        for( size_t i = 0; i < ss.functionIdMap().at( 1 ).size(); i++ )
        {
            geom.row( ss.functionIdMap().at( 1 ).at( i ) ) = greville_2.row( i ) + Eigen::RowVector2d( 1, 0 );
        }
        return geom;
    }();

    // Make an evaluator to test the extraction operators
    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto test_center_pt = [&]( const Face& v, const Eigen::Vector2d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas.parentDomain( v ), parentDomainBoundary( atlas, v ) );
        evaler.localizeElement( v );
        evaler.localizePoint( ppt );

        const Eigen::Vector2d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Face( 0 ), Eigen::Vector2d( 0.5, 0.25 ) );
    test_center_pt( Face( 4 ), Eigen::Vector2d( 0.5, 0.75 ) );
    test_center_pt( Face( 40 ), Eigen::Vector2d( 1.5, 0.25 ) );
    test_center_pt( Face( 64 ), Eigen::Vector2d( 1.25, 0.625 ) );
    test_center_pt( Face( 68 ), Eigen::Vector2d( 1.75, 0.625 ) );
    test_center_pt( Face( 72 ), Eigen::Vector2d( 1.25, 0.875 ) );
    test_center_pt( Face( 76 ), Eigen::Vector2d( 1.75, 0.875 ) );
}

TEST_CASE( "Simple 3d multi patch hierarchical spline space" )
{
    const KnotVector kv1( {0, 0, 0, 1, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0, 0.5, 1, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1}, 1e-10 );

    const size_t degree = 2;
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv1, kv1, kv2 }, { degree, degree, degree } ) );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv2, kv2, kv4 }, { degree, degree, degree } ) );

    const auto mp_ss_1 = std::make_shared<const MultiPatchSplineSpace>( buildMultiPatchSplineSpace( { tp_ss_1, tp_ss_1 }, {
        { { 0, Dart( 7 ) }, { 1, Dart( 19 ) } }
    } ) );

    const auto mp_ss_2 = std::make_shared<const MultiPatchSplineSpace>( buildMultiPatchSplineSpace( { tp_ss_2, tp_ss_2 }, {
        { { 0, Dart( 7 ) }, { 1, Dart( 19 ) } }
    } ) );

    const HierarchicalMultiPatchSplineSpace ss = buildHierarchicalSplineSpace( { mp_ss_1, mp_ss_2 }, {
        { Volume( 24 ), Volume( 48 ), Volume( 72 ) },
        { Volume( 0 ), Volume( 24 ), Volume( 48 ), Volume( 72 ), Volume( 96 ), Volume( 120 ), Volume( 144 ), Volume( 168 ) }
    } );

    CHECK( cellCount( ss.basisComplex().parametricAtlas().cmap(), 3 ) == 11 );
    CHECK( cellCount( ss.basisComplex().parametricAtlas().cmap(), 2 ) == 50 );
    CHECK( cellCount( ss.basisComplex().parametricAtlas().cmap(), 1 ) == 75 );
    CHECK( cellCount( ss.basisComplex().parametricAtlas().cmap(), 0 ) == 37 );
    CHECK( ss.numFunctions() == 78 );

    iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const Volume& v ) {
        const auto [level, dart] = ss.basisComplex().parametricAtlas().cmap().unrefinedAncestorDartOfCell( v );
        if( level == 0 ) CHECK( ss.connectivity( v ).size() == 27 );
        else CHECK( ss.connectivity( v ).size() >= 27 );
        CHECK( util::equals( ss.extractionOperator( v ).colwise().sum().transpose(), Eigen::VectorXd::Ones( 27 ), 1e-9 ) );
        return true;
    } );

    const Eigen::MatrixX3d geom = [&]() {
        const Eigen::MatrixXd greville_1 = grevillePoints( *ss.constituents().at( 0 ) );
        const Eigen::MatrixXd greville_2 = grevillePoints( *ss.constituents().at( 1 ) );

        Eigen::MatrixXd geom = Eigen::MatrixXd::Zero( ss.numFunctions(), 3 );

        for( size_t i = 0; i < ss.functionIdMap().at( 0 ).size(); i++ )
        {
            geom.row( ss.functionIdMap().at( 0 ).at( i ) ) = greville_1.row( i );
        }
        for( size_t i = 0; i < ss.functionIdMap().at( 1 ).size(); i++ )
        {
            geom.row( ss.functionIdMap().at( 1 ).at( i ) ) = greville_2.row( i ) + Eigen::RowVector2d( 1, 0 );
        }
        return geom;
    }();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_3d_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto test_center_pt = [&]( const Volume& v, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas.parentDomain( v ), parentDomainBoundary( atlas, v ) );
        evaler.localizeElement( v );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Volume( 26 ), Eigen::Vector3d( 0.5, 0.5, 0.75 ) );
    test_center_pt( Volume( 48 ), Eigen::Vector3d( 0.25, 0.25, 0.125 ) );
    test_center_pt( Volume( 72 ), Eigen::Vector3d( 0.75, 0.25, 0.125 ) );
    test_center_pt( Volume( 96 ), Eigen::Vector3d( 0.25, 0.75, 0.125 ) );
    test_center_pt( Volume( 120 ), Eigen::Vector3d( 0.75, 0.75, 0.125 ) );
    test_center_pt( Volume( 144 ), Eigen::Vector3d( 0.25, 0.25, 0.375 ) );
    test_center_pt( Volume( 168 ), Eigen::Vector3d( 0.75, 0.25, 0.375 ) );
    test_center_pt( Volume( 192 ), Eigen::Vector3d( 0.25, 0.75, 0.375 ) );
    test_center_pt( Volume( 216 ), Eigen::Vector3d( 0.75, 0.75, 0.375 ) );
    test_center_pt( Volume( 432 ), Eigen::Vector3d( 1.5, 0.5, 0.25 ) );
    test_center_pt( Volume( 456 ), Eigen::Vector3d( 1.5, 0.5, 0.75 ) );

}