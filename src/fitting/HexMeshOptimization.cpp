#include <HexMeshOptimization.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <HexOpt.h>

namespace fitting
{

using namespace HexOpt;
VertexPositionsFunc optimizeMesh( const topology::TetMeshCombinatorialMap& tet_mesh,
                                  const VertexPositionsFunc& tet_positions,
                                  const topology::CombinatorialMap& hex_mesh,
                                  const VertexPositionsFunc& hex_positions,
                                  const std::string& outputfile,
                                  const bool move_bdry_points )
{
    using namespace topology;
    std::map<Dart::IndexType, size_t> hex_vert_ids;
    HexMeshInfo hexMeshInfo = [&](){
        std::vector<std::vector<int>> hexes;
        std::vector<double> flattened_hex_positions;
        std::vector<std::vector<int>> adjacency_list( cellCount( hex_mesh, 0 ) );

        size_t i = 0;
        iterateCellsWhile( hex_mesh, 0, [&]( const Vertex& v ) {
            const Eigen::Vector3d pos = hex_positions( v );
            flattened_hex_positions.push_back( pos[0] );
            flattened_hex_positions.push_back( pos[1] );
            flattened_hex_positions.push_back( pos[2] );
            hex_vert_ids.emplace( lowestDartId( hex_mesh, v ), i++ );
            return true;
        } );
        // std::cout << flattened_hex_positions << std::endl;

        const IndexingFunc hex_indexing = [&]( const Vertex& v ) {
            return hex_vert_ids.at( lowestDartId( hex_mesh, v ) );
        };
        iterateCellsWhile( hex_mesh, 3, [&]( const Cell& c ) {
            hexes.emplace_back();
            std::vector<int>& hex = hexes.back();
            hex.reserve( 8 );
            hex.push_back( 3*hex_indexing( Vertex( c.dart() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, 1, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {1, 1}, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {1, 1, 1}, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {2, 1, 1, 2, 1}, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {2, 1, 1, 2}, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {2, 1, 1, 2, -1}, c.dart() ).value() ) ) );
            hex.push_back( 3*hex_indexing( Vertex( phi( hex_mesh, {2, 1, 1, 2, -1, -1}, c.dart() ).value() ) ) );
            // std::cout << "Hex: " << hex << std::endl;
            return true;
        } );

        iterateCellsWhile( hex_mesh, 1, [&]( const Edge& e ) {
            const Vertex v1( e.dart() );
            const Vertex v2( phi( hex_mesh, 1, e.dart() ).value() );
            const int idx1 = hex_indexing( v1 );
            const int idx2 = hex_indexing( v2 );
            adjacency_list[idx1].push_back( idx2 );
            adjacency_list[idx2].push_back( idx1 );
            return true;
        } );

        {
            CombinatorialMapBoundary bdry( hex_mesh );
            std::cout << "Boundary vertices: " << cellCount( bdry, 0 ) << std::endl;
        }

        return HexMeshInfo{ static_cast<int>( cellCount( hex_mesh, 0 ) ),
                            static_cast<int>( cellCount( hex_mesh, 3 ) ),
                            hexes,
                            flattened_hex_positions,
                            adjacency_list };
    }();
    FeatureInfo featureInfo{ hexMeshInfo.pNum };

    const TriMeshInfo triMeshInfo = [&](){
        CombinatorialMapBoundary bdry( tet_mesh );
        const auto tri_positions = boundaryVertexPositions( bdry, tet_positions );

        std::vector<std::vector<int>> tris;
        std::vector<std::vector<double>> tri_points;
        std::vector<std::vector<std::vector<double>>> triEdgeX;

        size_t i = 0;
        std::map<Dart::IndexType, size_t> tri_vert_ids;
        iterateCellsWhile( bdry, 0, [&]( const Vertex& v ) {
            const Eigen::Vector3d pos = tri_positions( v );
            tri_points.push_back( { pos[0], pos[1], pos[2] } );
            tri_vert_ids.emplace( lowestDartId( bdry, v ), i++ );
            return true;
        } );

        const IndexingFunc tri_indexing = [&]( const Vertex& v ) {
            return tri_vert_ids.at( lowestDartId( bdry, v ) );
        };
        iterateCellsWhile( bdry, 2, [&]( const Face& c ) {
            tris.emplace_back();
            std::vector<int>& tri = tris.back();
            tri.reserve( 8 );
            tri.push_back( tri_indexing( Vertex( c.dart() ) ) );
            tri.push_back( tri_indexing( Vertex( phi( bdry, 1, c.dart() ).value() ) ) );
            tri.push_back( tri_indexing( Vertex( phi( bdry, {1, 1}, c.dart() ).value() ) ) );

            triEdgeX.push_back({ { tri_points[tri[1]][0] - tri_points[tri[0]][0], tri_points[tri[1]][1] - tri_points[tri[0]][1], tri_points[tri[1]][2] - tri_points[tri[0]][2] },
                                 { tri_points[tri[2]][0] - tri_points[tri[0]][0], tri_points[tri[2]][1] - tri_points[tri[0]][1], tri_points[tri[2]][2] - tri_points[tri[0]][2] },
                                 { tri_points[tri[2]][0] - tri_points[tri[1]][0], tri_points[tri[2]][1] - tri_points[tri[1]][1], tri_points[tri[2]][2] - tri_points[tri[1]][2] } });
            return true;
        } );
        return TriMeshInfo{
            static_cast<int>( cellCount( bdry, 2 ) ), tris, tri_points, triEdgeX };
    }();

    optimizeMesh( hexMeshInfo, featureInfo, triMeshInfo, outputfile, move_bdry_points );

    return [&hex_mesh, hex_vert_ids, flat_points = hexMeshInfo.x]( const Vertex& v ) -> Eigen::Vector3d {
        auto it = hex_vert_ids.find( lowestDartId( hex_mesh, v ) );
        if ( it != hex_vert_ids.end() ) {
            size_t idx = it->second;
            return Eigen::Vector3d( flat_points[3 * idx], flat_points[3 * idx + 1], flat_points[3 * idx + 2] );
        }
        throw std::runtime_error( "Vertex not found in optimized mesh positions." );
    };
}
}