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

        // Need to make a map that goes from vertex index to vertex position, and then
        // iterate over all vertices to pull out the positions in correct order
        std::map<int,const topology::Vertex*> vert_idx_to_vert;
        std::map<int,const Eigen::Vector3d> vert_to_position;
        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );

        std::vector<int> vert_idxs;

//        std::cout << " Got here 0" << std::endl;
        int max_vert = 0;
        iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& v ) {
            const int idx = bdry_vertex_ids(v);
            vert_idx_to_vert.insert({idx, &v});
            if ( idx > max_vert )
            {
                max_vert = idx;
            }
            vert_idxs.push_back(idx);
            const Eigen::Vector3d point = bdry_positions( v );
            vert_to_position.insert({idx, point});
            return true;
        } );

//        std::cout << vert_idxs << std::endl;
//        std::sort (vert_idxs.begin(), vert_idxs.end());
//        std::cout << vert_idxs << std::endl;

//        std::cout << " Got here 1" << std::endl;
        for (int i = 0; i < max_vert+1; ++i)
        {
//            std::cout << "a" << std::endl;
            const topology::Vertex* v = vert_idx_to_vert.at(i);
//            std::cout << "b";
            const topology::Vertex& v1 = *v;
//            std::cout << "c";
            const Eigen::Vector3d point = vert_to_position.at(i);
//            std::cout << "d";
            file << "v " << point( 0 ) << " " << point( 1 ) << " " << point( 2 ) << std::endl;
        }


//        std::cout << " Got here 2" << std::endl;

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
