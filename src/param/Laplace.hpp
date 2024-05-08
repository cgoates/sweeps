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

Eigen::VectorXd sweepEmbedding( const topology::TetMeshCombinatorialMap& map,
                                const std::vector<bool>& zero_bcs,
                                const std::vector<bool>& one_bcs,
                                const std::vector<Normal>& normals );

Eigen::VectorXd sweepEmbedding( const topology::CombinatorialMap& map,
                                const std::function<double( const topology::Edge& )>& edge_weights,
                                const std::vector<bool>& zero_bcs,
                                const std::vector<bool>& one_bcs );

Eigen::MatrixX2d tutteEmbedding( const topology::CombinatorialMap& map,
                                 const VertexPositionsFunc& vert_positions,
                                 const std::function<std::optional<Eigen::Vector2d>( const topology::Vertex& )>& constraints,
                                 const bool shape_preserving = true );

Eigen::MatrixXd solveLaplaceSparse(
    const topology::CombinatorialMap& map,
    const std::function<double( const topology::Edge& )>& edge_weights,
    const std::function<std::optional<Eigen::VectorXd>( const topology::Vertex& )>& constraints,
    const size_t n_constrained_verts,
    const size_t data_dim );