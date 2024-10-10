#pragma once
#include <string>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>

namespace io
{
    /// @brief write out the boundary of a tet mesh to OBJ file format
    /// @param to_output The tet mesh
    /// @param filename The output filename, including the extension.
    void outputTetMeshBoundaryToOBJ( const topology::TetMeshCombinatorialMap& to_output, const std::string& filename );
} // namespace io
