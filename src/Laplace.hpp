#pragma once
#include <Eigen/Dense>
#include <set>

class VertexId;
class Normal;

namespace topology
{
    class TetMeshCombinatorialMap;
}

Eigen::VectorXd solveLaplaceSparse( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals );