#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriMeshCombinatorialMap.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentDomain.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>

TEST_CASE( "Simple circle evaluations" )
{
    // This but fit boundary edges to a circle:
    //   *
    //  /|\.
    // *-*-*
    //  \|/
    //   *
    const SimplicialComplex mesh{
        { { 0, 1, 4 }, { 0, 4, 2 }, { 4, 1, 3 }, { 4, 3, 2 } },
        { { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.2 }, { 0.0, 0.0, 0.0 } } };

    const topology::TriMeshCombinatorialMap cmap( mesh );
    const auto vert_id = indexingOrError( cmap, 0 );
    const param::TriangleParametricAtlas param( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            if( boundaryAdjacent( cmap, e ) )
                CHECK( util::equals( pos.norm(), 1.0, 1e-12 ) );
            else
                CHECK( util::equals( pos.norm(), 0.5, 1e-12 ) );
            return true;
        } );

        iterateAdjacentCells( cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param.parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const Eigen::Vector2d expected = vert_positions( v );
            CHECK( util::equals( pos, expected, 1e-12 ) );
            return true;
        } );

        const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( param, f ) );
        const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
        CHECK( util::equals( pos.norm(), 2.0/3.0, 1e-12 ) );
        return true;
    } );
}

TEST_CASE( "Circle evaluation with non-boundary-adjacent cells" )
{
    // An equilateral triangle within a circle
    std::vector<Eigen::Vector3d> points;
    for( const double angle : { 0.0, 2.0 * std::numbers::pi / 3.0, 4.0 * std::numbers::pi / 3.0 } )
    {
        points.emplace_back( 0.5 * std::cos( angle ), 0.5 * std::sin( angle ), 0.0 );
        points.emplace_back( std::cos( angle ), std::sin( angle ), 0.0 );
    }
    const SimplicialComplex mesh{
        { {0, 2, 4}, {0, 1, 2}, {2, 1, 3}, {2, 3, 4}, {4, 3, 5}, {4, 5, 0}, {0, 5, 1} },
        points
    };

        const topology::TriMeshCombinatorialMap cmap( mesh );
    const auto vert_id = indexingOrError( cmap, 0 );
    const param::TriangleParametricAtlas param( cmap );
    const auto vert_positions = [&]( const topology::Vertex& v ) {
        return mesh.points.at( vert_id( v ) ).head<2>();
    };
    const mapping::TriangleMeshCircleMapping geom_mapping( param, vert_positions );

    const param::ParentDomain pd = param::simplexDomain( 2 );
    iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
        iterateAdjacentCells( cmap, f, 1, [&]( const topology::Edge& e ) {
            const param::ParentPoint ppt = pointOnBoundary( pd, parentDomainBoundary( param, e ) );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            if( boundaryAdjacent( cmap, e ) )
                CHECK( util::equals( pos.norm(), 1.0, 1e-12 ) );
            return true;
        } );

        iterateAdjacentCells( cmap, f, 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint ppt = param.parentPoint( v );
            const Eigen::Vector2d pos = geom_mapping.evaluate( f, ppt );
            const Eigen::Vector2d expected = vert_positions( v );
            CHECK( util::equals( pos, expected, 1e-12 ) );
            return true;
        } );
        return true;
    } );

    const param::ParentPoint ppt = pointOnBoundary( pd, {false, false, false} );
    const Eigen::Vector2d pos = geom_mapping.evaluate( topology::Face( topology::Dart( 0 ) ), ppt );
    CHECK( util::equals( pos, Eigen::Vector2d( 0.0, 0.0 ), 1e-12 ) );
}
