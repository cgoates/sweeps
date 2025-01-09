#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentDomain.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

TEST_CASE( "Closest point test" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  / \.
    // *---*
    //  \ /
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 3 }, { 0, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 } } };

    const auto cmap = std::make_shared<const topology::TriMeshCombinatorialMap>( mesh );
    const auto vert_id = indexingOrError( *cmap, 0 );
    const auto param = std::make_shared<const param::TriangleParametricAtlas>( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) );
    };
    const mapping::TriangleMeshMapping geom_mapping( param, vert_positions, 3 );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( *cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( *cmap, f, 1, [&]( const topology::Edge& e ) {
            if( not boundaryAdjacent( *cmap, e ) ) return true;
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( *param, e ) );
            const Eigen::Vector3d pos = geom_mapping.evaluate( f, ppt );
            const auto ppt_rt = geom_mapping.closestPoint( pos + Eigen::Vector3d( 0, 0, 1 ) );
            CHECK( ppt_rt.second.mBaryCoordIsZero == ppt.mBaryCoordIsZero );
            CHECK( util::equals( ppt_rt.second.mPoint, ppt.mPoint, 1e-10 ) );
            CHECK( f == ppt_rt.first );
            return true;
        } );
        return true;
    } );

    const auto test_point = [&]( const Eigen::Vector3d& point, const topology::Face& expected_f, const param::ParentPoint& expected_ppt ) {
        const auto ppt_rt = geom_mapping.closestPoint( point );

        CHECK( ppt_rt.second.mBaryCoordIsZero == expected_ppt.mBaryCoordIsZero );
        CHECK( util::equals( ppt_rt.second.mPoint, expected_ppt.mPoint, 1e-10 ) );
        CHECK( expected_f == ppt_rt.first );
    };

    const auto test_point_roundtrip = [&]( const Eigen::Vector3d& point, const Eigen::Vector3d& expected_point ) {
        const auto ppt_rt = geom_mapping.closestPoint( point );
        const auto point_rt = geom_mapping.evaluate( ppt_rt.first, ppt_rt.second );

        CHECK( util::equals( expected_point, point_rt, 1e-10 ) );
    };

    test_point( { 0.0, 0.1, 0.1 }, topology::Face( 3 ), param::ParentPoint( pd, Eigen::Vector2d( 0.45, 0.1 ), {false, false, false} ) );
    test_point( { 0.0, 5.0, 0.1 }, topology::Face( 3 ), param::ParentPoint( pd, Eigen::Vector2d( 0.0, 1 ), {true, true, false} ) );
    test_point_roundtrip( { 10.0, 0.01, 0.0 }, { 1.0, 0.0, 0.0 } );
    test_point_roundtrip( { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 } );
    test_point( { 1.0, -1.0, 0.0 }, topology::Face( 0 ), param::ParentPoint( pd, Eigen::Vector2d( 0.5, 0.5 ), {true, false, false} ) );
}