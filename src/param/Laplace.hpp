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
    class Edge;
}

Eigen::VectorXd solveLaplaceSparse( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals );

Eigen::VectorXd solveLaplaceSparse( const topology::CombinatorialMap& map,
                                    const std::function<double(const topology::Edge&)>& edge_weights,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs );

Eigen::VectorXd solveLaplaceSparse( const topology::CombinatorialMap& map,
                                    const std::function<double( const topology::Edge& )>& edge_weights,
                                    const std::function<std::optional<double>( const topology::Vertex& )>& constraints,
                                    const size_t n_constrained_verts );