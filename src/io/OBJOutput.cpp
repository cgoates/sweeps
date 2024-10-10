#include <OBJOutput.hpp>
#include <fstream>
#include <sstream>
#include <Simplex.hpp>
#include <chrono>
#include <ctime>

namespace io
{
    void outputTetMeshBoundaryToOBJ( const topology::TetMeshCombinatorialMap& to_output, const std::string& filename )
    {
        const auto vertex_ids = indexingOrError( to_output, 0 );
        const auto vertex_positions = [&]( const topology::Vertex& v ) {
            return to_output.simplicialComplex().points.at( vertex_ids( v ) );
        };

        // extract the boundary triangular mesh
        const topology::CombinatorialMapBoundary bdry( to_output );

        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
        const auto bdry_positions = [&]( const topology::Vertex& v ) {
            return vertex_positions( bdry.toUnderlyingCell( v ) );
        };

        if( cellCount( bdry, 2 ) == 0 )
        {
            throw std::invalid_argument( "A combinatorial map with no boundary was given for boundary OBJ output." );
        }

        std::ofstream file;
        file.open( filename );

        const auto start = std::chrono::system_clock::now();
        const std::time_t start_time = std::chrono::system_clock::to_time_t( start );

        file << "# File created on " << std::ctime( &start_time );

        // Map from cmap vertex index to obj vertex index, which is the position in the OBJ file vertex list
        std::map<size_t, size_t> cmap_vids_to_obj_vids;

        size_t vert_ii = 1;// OBJ files use one-based indexing
        iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& v ) {
            const Eigen::Vector3d point = bdry_positions( v );
            file << "v " << point.transpose() << std::endl;
            cmap_vids_to_obj_vids.insert( { bdry_vertex_ids( v ), vert_ii++ } );
            return true;
        } );

        file << std::endl;

        iterateCellsWhile( bdry, 2, [&]( const topology::Face& f ) {
            file << "f";
            topology::Dart current_d = f.dart();
            do
            {
                const topology::Vertex v( current_d );
                file << " " << cmap_vids_to_obj_vids.at( bdry_vertex_ids( v ) );
                current_d = phi( bdry, 1, current_d ).value();
            } while( current_d != f.dart() );
            file << std::endl;

            return true;
        } );

        file.close();
    }

} // namespace io
