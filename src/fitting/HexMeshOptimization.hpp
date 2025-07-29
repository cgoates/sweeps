#pragma once
#include <VertexPositionsFunc.hpp>
#include <TetMeshCombinatorialMap.hpp>

namespace fitting
{
    VertexPositionsFunc optimizeMesh( const topology::TetMeshCombinatorialMap& tet_mesh,
                                      const VertexPositionsFunc& tet_positions,
                                      const topology::CombinatorialMap& hex_mesh,
                                      const VertexPositionsFunc& hex_positions,
                                      const std::string& outputfile,
                                      const bool move_bdry_points = false );
}