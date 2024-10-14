#include <catch2/catch_test_macros.hpp>
#include <Dijkstra.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <VTKOutput.hpp>
#include <SimplexUtilities.hpp>

using namespace topology;

TEST_CASE( "Shortest path finding" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/dodec.inp", "Surface1", "Surface2" );
    const topology::TetMeshCombinatorialMap cmap( sweep_input.mesh );
    const topology::CombinatorialMapBoundary bdry( cmap );

    const auto vert_ids = indexingOrError( bdry, 0 );
    const auto bdry_positions = [&]( const topology::Vertex& v ){
        return sweep_input.mesh.points.at( vert_ids( v ) );
    };

    topology::Vertex start_v( 0 );
    topology::Vertex end_v( 0 );

    const auto test_equal_vertices = [&vert_ids]( const topology::Vertex& end_v ) {
        return [&]( const topology::Vertex& test_v ){
            return vert_ids( test_v ) == vert_ids( end_v );
        };
    };

    const auto path_length = [&bdry,&bdry_positions]( const std::vector<topology::Edge>& path ) {
        return std::transform_reduce( path.begin(), path.end(), 0.0, std::plus<>(), [&]( const topology::Edge& e ){ return edgeLength( bdry, bdry_positions, e ); } );
    };

    const auto nearby_vertices = [&]( const std::vector<topology::Edge>& edges ) {
        std::set<topology::Vertex> out;
        for( const auto& edge : edges )
        {
            iterateDartsOfCell( bdry, edge, [&]( const topology::Face& f ) {
                iterateAdjacentCells( bdry, f, 0, [&]( const topology::Vertex& v ) {
                    out.insert( v );
                    return true;
                } );
                return true;
            } );
        }
        return out;
    };

    iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& start_v ) {
        iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& end_v ) {
            if( test_equal_vertices( start_v )( end_v ) ) return true;

            const auto edges = topology::shortestPath( bdry, bdry_positions, start_v, test_equal_vertices( end_v ) );
            const double shortest_path_length = path_length( edges );

            for( const topology::Vertex& mid_v : nearby_vertices( edges ) )
            {
                if( test_equal_vertices( start_v )( mid_v ) ) continue;
                if( test_equal_vertices( mid_v )( end_v ) ) continue;
                const auto to_mid_edges = topology::shortestPath( bdry, bdry_positions, start_v, test_equal_vertices( mid_v ) );
                const double to_mid_path_length = path_length( to_mid_edges );
                const auto from_mid_edges = topology::shortestPath( bdry, bdry_positions, mid_v, test_equal_vertices( end_v ) );
                const double from_mid_path_length = path_length( from_mid_edges );

                CHECK( shortest_path_length < to_mid_path_length + from_mid_path_length + 1e-9 );
            }
            return true;
        } );
        return true;
    } );
}