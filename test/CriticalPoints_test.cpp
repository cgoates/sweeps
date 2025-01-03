#include <catch2/catch_test_macros.hpp>
#include <SweepInput.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Logging.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>
#include <CriticalPoints.hpp>

const auto test_critical_points = []( const Eigen::VectorXd& vals, const int expected_euler ) {
    const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();
    
    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );

    const auto vert_ids = indexingOrError( map, 0 );

    CHECK( expected_euler == reparam::lowerLinkEulerCharacteristic( map, map.vertexOfId( 8 ), [&]( const topology::Vertex& v ) { return vals( vert_ids( v ) ); } ) );
    if( expected_euler == 1 )
        CHECK( not reparam::isCriticalPoint( map, map.vertexOfId( 8 ), [&]( const topology::Vertex& v ) { return vals( vert_ids( v ) ); } ) );
    else
        CHECK( reparam::isCriticalPoint( map, map.vertexOfId( 8 ), [&]( const topology::Vertex& v ) { return vals( vert_ids( v ) ); } ) );
};

TEST_CASE( "Critical points test", "" )
{
    test_critical_points( ( Eigen::VectorXd(12) << 0, 0, 0, 0, 1, 1, 1, 1, 0.5 ).finished(), 1 );
    test_critical_points( ( Eigen::VectorXd(12) << 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.25 ).finished(), 2 );
    test_critical_points( ( Eigen::VectorXd(12) << 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.75 ).finished(), 0 );
    test_critical_points( ( Eigen::VectorXd(12) << 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.75 ).finished(), 2 );
    test_critical_points( ( Eigen::VectorXd(12) << 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.25 ).finished(), 0 );
}
