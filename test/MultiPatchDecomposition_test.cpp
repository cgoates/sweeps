#include <catch2/catch_test_macros.hpp>
#include <TPCombinatorialMap.hpp>
#include <MultiPatchDecomposition.hpp>
#include <QuadMeshCombinatorialMap.hpp>
#include <DelaunayTriangulation.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <MeshInput.hpp>
#include <VTKOutput.hpp>
#include <GlobalCellMarker.hpp>

using namespace topology;
TEST_CASE( "Testing corner identification for multipatch" )
{
    const auto [faces, points] = io::loadOBJFile( SRC_HOME "/test/data/capsule_quadmesh.obj" );
    std::vector<std::array<VertexId, 4>> quads;
    quads.reserve( faces.size() );
    for( const auto& f : faces )
    {
        if( f.size() == 4 )
        {
            quads.push_back( { VertexId( f.at( 0 ) ), VertexId( f.at( 1 ) ), VertexId( f.at( 2 ) ), VertexId( f.at( 3 ) ) } );
        }
        else
        {
            throw std::runtime_error( "Only quads are supported in this test." );
        }
    }
    const std::shared_ptr<const QuadMeshCombinatorialMap> cmap = std::make_shared<const QuadMeshCombinatorialMap>( quads, points.size() );

    const auto vids = indexingOrError( *cmap, 0 );

    const auto corners = multiPatchCorners( *cmap );
    const auto vert_positions = [&]( const Vertex& v ) -> Eigen::Vector3d {
        return points.at( vids( v ) );
    };

    SimplicialComplex corner_points;
    iterateCellsWhile( *cmap, 0, [&]( const Vertex& v ) {
        if( corners.isMarked( v ) )
        {
            corner_points.points.push_back( points.at( vids( v ) ) );
            corner_points.simplices.emplace_back( corner_points.points.size() - 1 );
        }
        return true;
    } );
    std::cout << "Number of corner points: " << corner_points.points.size() << std::endl;

    const DelaunayTriangulation triangulation( cmap, vert_positions );

    io::VTKOutputObject output( corner_points );
    io::outputSimplicialFieldToVTK( output, "corner_points.vtu" );

    io::outputCMap( triangulation, vert_positions, "triangulation.vtu" );

    SimplicialComplex edges;
    edges.points = points;
    iterateCellsWhile( *cmap, 1, [&]( const Edge& e ) {
        edges.simplices.emplace_back( vids( Vertex( e.dart() ) ),
                                      vids( Vertex( phi( *cmap, 1, e.dart() ).value() ) ) );
        return true;
    } );
    io::VTKOutputObject edge_output( edges );
    io::outputSimplicialFieldToVTK( edge_output, "edges.vtu" );

    const MultiPatchDecomposition mp_decomp = multiPatchDecomposition( *cmap );

    const MultiPatchCombinatorialMap mp_cmap( mp_decomp.constituents, mp_decomp.connections );

    CHECK( cellCount( mp_cmap, 2 ) == cellCount( *cmap, 2 ) );
    CHECK( cellCount( mp_cmap, 1 ) == cellCount( *cmap, 1 ) );
    CHECK( cellCount( mp_cmap, 0 ) == cellCount( *cmap, 0 ) );
}

TEST_CASE( "Testing multipatch decomposition 2" )
{
    const auto [faces, points] = io::loadOBJFile( SRC_HOME "/test/data/hookBase.obj" );
    std::vector<std::array<VertexId, 4>> quads;
    quads.reserve( faces.size() );
    for( const auto& f : faces )
    {
        if( f.size() == 4 )
        {
            quads.push_back( { VertexId( f.at( 0 ) ), VertexId( f.at( 1 ) ), VertexId( f.at( 2 ) ), VertexId( f.at( 3 ) ) } );
        }
        else
        {
            throw std::runtime_error( "Only quads are supported in this test." );
        }
    }
    const std::shared_ptr<const QuadMeshCombinatorialMap> cmap = std::make_shared<const QuadMeshCombinatorialMap>( quads, points.size() );

    const auto vids = indexingOrError( *cmap, 0 );

    const auto corners = multiPatchCorners( *cmap );
    const auto vert_positions = [&]( const Vertex& v ) -> Eigen::Vector3d {
        return points.at( vids( v ) );
    };

    SimplicialComplex corner_points;
    iterateCellsWhile( *cmap, 0, [&]( const Vertex& v ) {
        if( corners.isMarked( v ) )
        {
            corner_points.points.push_back( points.at( vids( v ) ) );
            corner_points.simplices.emplace_back( corner_points.points.size() - 1 );
        }
        return true;
    } );
    std::cout << "Number of corner points: " << corner_points.points.size() << std::endl;

    const DelaunayTriangulation triangulation( cmap, vert_positions );

    io::VTKOutputObject output( corner_points );
    io::outputSimplicialFieldToVTK( output, "corner_points.vtu" );

    io::outputCMap( triangulation, vert_positions, "triangulation.vtu" );

    SimplicialComplex edges;
    edges.points = points;
    iterateCellsWhile( *cmap, 1, [&]( const Edge& e ) {
        edges.simplices.emplace_back( vids( Vertex( e.dart() ) ),
                                      vids( Vertex( phi( *cmap, 1, e.dart() ).value() ) ) );
        return true;
    } );
    io::VTKOutputObject edge_output( edges );
    io::outputSimplicialFieldToVTK( edge_output, "edges.vtu" );

    const MultiPatchDecomposition mp_decomp = multiPatchDecomposition( *cmap );

    CHECK( mp_decomp.constituents.size() == 18 );

    const MultiPatchCombinatorialMap mp_cmap( mp_decomp.constituents, mp_decomp.connections );

    CHECK( cellCount( mp_cmap, 2 ) == cellCount( *cmap, 2 ) );
    CHECK( cellCount( mp_cmap, 1 ) == cellCount( *cmap, 1 ) );
    CHECK( cellCount( mp_cmap, 0 ) == cellCount( *cmap, 0 ) );
}


