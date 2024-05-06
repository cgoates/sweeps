#pragma once
#include <Eigen/Dense>
#include <set>
#include <VertexPositionsFunc.hpp>

class VertexId;
class Normal;

namespace topology
{
    class TetMeshCombinatorialMap;
    class CombinatorialMap;
}

Eigen::VectorXd solveLaplaceSparse( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals );

Eigen::VectorXd solveLaplaceSparse( const topology::CombinatorialMap& map,
                                    const VertexPositionsFunc& v_positions,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals );