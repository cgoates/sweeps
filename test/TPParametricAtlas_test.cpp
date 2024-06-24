#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <TPParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

using namespace topology;
using namespace param;
TEST_CASE( "3x4 TP combinatorial map" )
{
    const CombinatorialMap1d source_topo( 3 );
    const CombinatorialMap1d line_topo( 4 );
    const TPCombinatorialMap tp_topo( source_topo, line_topo );

    const ParametricAtlas1d source( source_topo, Eigen::Vector3d( 1, 2, 3 ) );
    const ParametricAtlas1d line( line_topo, Eigen::Vector4d( 0.1, 0.2, 0.3, 0.4 ) );
    const TPParametricAtlas tp( tp_topo, source, line );

    const ParentDomain expected_domain = cubeDomain( 2 );
    iterateCellsWhile( tp_topo, 2, [&]( const Cell& c ) {
        CHECK( tp.parentDomain( c ) == expected_domain );
        return true;
    } );

    std::array<ParentPoint, 4> expected_points( {
        ParentPoint( expected_domain, Eigen::Vector2d( 0.0, 0.0 ), { false, true, false, true } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 1.0, 0.0 ), { true, false, false, true } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 1.0, 1.0 ), { true, false, true, false } ),
        ParentPoint( expected_domain, Eigen::Vector2d( 0.0, 1.0 ), { false, true, true, false } )
    } );
    iterateDartsWhile( tp_topo, [&]( const Dart& d ) {
        const auto ppt = tp.parentPoint( d );
        const auto& expected_ppt = expected_points.at( d.id() % 4 );
        CHECK( ppt.mDomain == expected_domain );
        CHECK( ppt.mBaryCoordIsZero == expected_ppt.mBaryCoordIsZero );
        CHECK( util::equals( ppt.mPoint, expected_ppt.mPoint, 1e-10 ) );
        return true;
    } );

    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 0 ), 2 ) ), Eigen::Vector2d( 1.0, 0.1 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 20 ), 2 ) ), Eigen::Vector2d( 3.0, 0.2 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 37 ), 2 ) ), Eigen::Vector2d( 1.0, 0.4 ), 1e-10 ) );
    CHECK( util::equals( tp.parametricLengths( Cell( Dart( 30 ), 2 ) ), Eigen::Vector2d( 2.0, 0.3 ), 1e-10 ) );
}
