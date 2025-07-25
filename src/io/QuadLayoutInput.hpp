#pragma once
#include <string>
#include <SweepInput.hpp>
#include <map>

namespace io
{
    /// \brief Rewrites the quad layout file in terms of vertex ids in the 3D mesh.
    /// Searches for a vertex in the 3D mesh that matches each vertex in the input OBJ file,
    /// then replaces the vertex ids in the barycentric file.
    void rewriteBaryCoordFile( const SimplicialComplex& mesh3d,
                               const std::string& obj_filename,
                               const std::string& bary_filename_in,
                               const std::string& bary_filename_out );

    /// \brief Loads curves consisting of barycentric coordinate in from a barycentric file.
    std::vector<std::vector<BarycentricPoint>> loadBaryCoords( const std::string& bary_filename_in );
} // namespace io