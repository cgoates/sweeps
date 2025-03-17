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
    /// @returns A map from the vertex index in the combinatorial map to the vertex index in the OBJ file (which is 1-based).
    std::map<size_t, size_t> outputTetMeshBoundaryToOBJ( const topology::TetMeshCombinatorialMap& to_output, const std::string& filename );
} // namespace io
