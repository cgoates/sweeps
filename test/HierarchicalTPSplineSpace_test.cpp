#include <catch2/catch_test_macros.hpp>
#include <HierarchicalTPSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <VTKOutput.hpp>
#include <SplineSpaceEvaluator.hpp>

using namespace topology;
using namespace param;
using namespace basis;

constexpr bool VTK_OUTPUT = false;

TEST_CASE( "Simple hierarchical atlas 1" )
{
    const KnotVector kv1( {0, 0, 0, 1, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0, 0.5, 1, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1}, 1e-10 );

    const size_t degree = 2;
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv1, kv2 }, { degree, degree } ) );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv2, kv4 }, { degree, degree } ) );

    const HierarchicalTPSplineSpace ss = buildHierarchicalSplineSpace( { tp_ss_1, tp_ss_2 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ), Face( Dart( 12 ) ) }
    } );

    const Eigen::MatrixX3d geom = ( Eigen::MatrixX3d( 17, 3 ) <<
        0, 0.5, 0,
        0.5, 0.5, 0,
        1, 0.5, 0,
        0, 1.5, 0,
        0.5, 1.5, 0,
        1, 1.5, 0,
        0, 2, 0,
        0.5, 2, 0,
        1, 2, 0,
        0, 0, 0,
        0.25, 0, 0,
        0.75, 0, 0,
        1, 0, 0,
        0, 0.25, 0,
        0.25, 0.25, 0,
        0.75, 0.25, 0,
        1, 0.25, 0 ).finished();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto test_center_pt = [&]( const Face& f, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas.parentDomain( f ), parentDomainBoundary( atlas, f ) );
        evaler.localizeElement( f );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Face( 5 ), Eigen::Vector3d( 0.5, 1.5, 0.0 ) );
    test_center_pt( Face( 8 ), Eigen::Vector3d( 0.25, 0.25, 0.0 ) );
    test_center_pt( Face( 12 ), Eigen::Vector3d( 0.75, 0.25, 0.0 ) );
}

TEST_CASE( "Simple hierarchical atlas 2" )
{
    const KnotVector kv1( {0, 0, 1, 1}, 1e-10 );
    const KnotVector kv2( {0, 0, 0.5, 1, 1}, 1e-10 );
    const KnotVector kv4( {0, 0, 0.25, 0.5, 0.75, 1, 1}, 1e-10 );
    const KnotVector kv8( { 0, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1 }, 1e-10 );

    const size_t degree = 1;
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv1, kv2 }, { degree, degree } ) );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv2, kv4 }, { degree, degree } ) );
    const auto tp_ss_3 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv4, kv8 }, { degree, degree } ) );

    const HierarchicalTPSplineSpace ss = buildHierarchicalSplineSpace( { tp_ss_1, tp_ss_2, tp_ss_3 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ) },
        { Face( Dart( 40 ) ), Face( Dart( 44 ) ), Face( Dart( 56 ) ), Face( Dart( 60 ) ) }
    } );

    const Eigen::MatrixX3d geom = ( Eigen::MatrixX3d( 12, 3 ) <<
        0, 1, 0,
        1, 1, 0,
        0, 2, 0,
        1, 2, 0,
        0, 0, 0,
        0.5, 0, 0,
        1, 0, 0,
        0, 0.5, 0,
        0.5, 0.5, 0,
        1, 0.5, 0,
        0.75, 0.75, 0,
        1.0, 0.75, 0 ).finished();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto test_center_pt = [&]( const Face& f, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas.parentDomain( f ), parentDomainBoundary( atlas, f ) );
        evaler.localizeElement( f );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );

        CHECK( util::equals( pt, expected, 1e-9 ) );
    };

    test_center_pt( Face( 5 ), Eigen::Vector3d( 0.5, 1.5, 0.0 ) );
    test_center_pt( Face( 8 ), Eigen::Vector3d( 0.25, 0.25, 0.0 ) );
    test_center_pt( Face( 12 ), Eigen::Vector3d( 0.75, 0.25, 0.0 ) );
    test_center_pt( Face( 80 ), Eigen::Vector3d( 0.625, 0.625, 0.0 ) );
}

TEST_CASE( "Simple hierarchical atlas 3" )
{
    const KnotVector kv1( { 0, 0, 1, 1}, 1e-10 );
    const KnotVector kv2( { 0, 0, 1, 2, 2}, 1e-10 );
    const KnotVector kv3 = nAdicRefine( kv1, 3 );
    const KnotVector kv6 = nAdicRefine( kv2, 3 );

    const size_t degree = 1;
    const auto tp_ss_1 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv1, kv2 }, { degree, degree } ) );
    const auto tp_ss_2 = std::make_shared<const TPSplineSpace>( buildBSpline( { kv3, kv6 }, { degree, degree } ) );

    const HierarchicalTPSplineSpace ss = buildHierarchicalSplineSpace( { tp_ss_1, tp_ss_2 }, {
        { Face( Dart( 4 ) ) },
        { Face( Dart( 0 ) ), Face( Dart( 4 ) ), Face( Dart( 8 ) ), Face( Dart( 12 ) ), Face( 16 ), Face( 20 ), Face( 24 ), Face( 28 ), Face( 32 ) }
    } );

    const Eigen::MatrixX3d geom = ( Eigen::MatrixX3d( 16, 3 ) <<
        0, 1, 0,
        1, 1, 0,
        0, 2, 0,
        1, 2, 0,
        0, 0, 0,
        0.3333, 0, 0,
        0.6666, 0, 0,
        1, 0, 0,
        0, 0.3333, 0,
        0.3333, 0.3333, 0,
        0.6666, 0.3333, 0,
        1, 0.3333, 0,
        0, 0.6666, 0,
        0.3333, 0.6666, 0,
        0.6666, 0.6666, 0,
        1, 0.6666, 0 ).finished();

    if( VTK_OUTPUT ) io::outputBezierMeshToVTK( ss, geom, "hier_bez_test.vtu" );

    eval::SplineSpaceEvaluator evaler( ss, 0 );

    const auto& atlas = ss.basisComplex().parametricAtlas();
    const auto test_center_pt = [&]( const Face& f, const Eigen::Vector3d& expected ) {
        const ParentPoint ppt = pointOnBoundary( atlas.parentDomain( f ), parentDomainBoundary( atlas, f ) );
        evaler.localizeElement( f );
        evaler.localizePoint( ppt );

        const Eigen::Vector3d pt = evaler.evaluateManifold( geom.transpose() );
        std::cout << pt.transpose() << std::endl;

        CHECK( util::equals( pt, expected, 1e-4 ) );
    };

    test_center_pt( Face( 5 ), Eigen::Vector3d( 0.5, 1.5, 0.0 ) );
    test_center_pt( Face( 8 ), Eigen::Vector3d( 1.0/6.0, 1.0/6.0, 0.0 ) );
    test_center_pt( Face( 12 ), Eigen::Vector3d( 0.5, 1.0/6.0, 0.0 ) );
    test_center_pt( Face( 32 ), Eigen::Vector3d( 1.0/6.0, 5.0/6.0, 0.0 ) );

    iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const Face f ){
        std::cout << f << std::endl << ss.connectivity( f ) << std::endl << ss.extractionOperator( f ) << std::endl << std::endl;
        return true;
    } );
}