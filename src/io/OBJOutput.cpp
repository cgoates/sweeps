#include <OBJOutput.hpp>
#include <fstream>
#include <sstream>
#include <Simplex.hpp>
#include <chrono>
#include <ctime>

namespace io
{
    void outputSimplicialMeshToOBJ( const topology::TetMeshCombinatorialMap& to_output, const std::string& filename )
    {
        const auto vertex_ids = indexingOrError( to_output, 0 );
        const auto vertex_positions = [&]( const topology::Vertex& v ) {
            return to_output.simplicialComplex().points.at( vertex_ids( v ) );
        };

        // extract the bounary triangular mesh
        topology::CombinatorialMapBoundary bdry(to_output);

        const auto bdry_positions = [&]( const topology::Vertex& v) {
            return vertex_positions( bdry.toUnderlyingCell( v ) );
        };

        const size_t n_cells = cellCount(bdry, 2);
        
        // ignore the empty case
        if (n_cells == 0)
        {
            std::cout << "An combinatorial map with empty boundary was given. No OBJ output was generated." << std::endl;
            return;
        }

        // only operate on simplicial objects of dimension 2
        // this is currently unnecessary since a Tet mesh is given
        const size_t dim = to_output.dim();
        if ( dim > 3 )
        {
            std::cout << "OBJ file format does not operate on combinatorial maps of dimension greater than 2" << std::endl;
            return;
        }
        else if ( dim < 3 )
        {
            std::cout << "Output of OBJ files for low-dimensional cell complexes is not currently supported" << std::endl;
            return;
        }

        std::ofstream file;
        file.open( filename );

        auto start = std::chrono::system_clock::now();
        std::time_t start_time = std::chrono::system_clock::to_time_t(start);

        file << "# File created on " << std::ctime(&start_time);

        iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& v ) {
            const Eigen::Vector3d point = bdry_positions( v );  
            file << "v " << point( 0 ) << " " << point( 1 ) << " " << point( 2 ) << std::endl;
            return true;
        } );

        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );


        iterateCellsWhile( bdry, 2, [&]( const topology::Face& f) {
            topology::Dart current_d = f.dart();
            file << "f ";
            do
            {
                const topology::Vertex v( current_d );
                file << bdry_vertex_ids(v) << " ";
                current_d = phi( bdry, 1, current_d ).value(); // value is okay here because you know that all phi1s are defined in a 2d cmap
            } while( current_d != f.dart() );
            file << std::endl;

            return true;
        } );

        file.close();
    }

} // namespace io
