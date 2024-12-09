#pragma once
#include <Eigen/Dense>
#include <set>
#include <VertexPositionsFunc.hpp>
#include <optional>

class VertexId;
class Normal;

namespace topology
{
    class TetMeshCombinatorialMap;
    class CombinatorialMap;
    class Edge;
}

namespace reparam
{
    enum class LaplaceEdgeWeights
    {
        Cotangent,
        InverseLength,
        VoronoiDual,
        BarycentricDual,
        Uniform
    };

    Eigen::VectorXd sweepEmbedding( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals,
                                    const LaplaceEdgeWeights& edge_weights = LaplaceEdgeWeights::Cotangent );

    Eigen::MatrixX2d
        tutteEmbedding( const topology::CombinatorialMap& map,
                        const VertexPositionsFunc& vert_positions,
                        const std::function<std::optional<Eigen::Vector2d>( const topology::Vertex& )>& constraints,
                        const bool shape_preserving = true );

    /// Defaults to embedding to a unit circle
    /// FIXME: Needs a way to specify orientation of the circle, not just randomly pick
    Eigen::MatrixX2d tutteEmbedding( const topology::CombinatorialMap& map, const VertexPositionsFunc& vert_positions );

    Eigen::MatrixXd
        solveLaplaceSparse( const topology::CombinatorialMap& map,
                            const std::function<double( const topology::Edge& )>& edge_weights,
                            const std::function<std::optional<Eigen::VectorXd>( const topology::Vertex& )>& constraints,
                            const size_t n_constrained_verts,
                            const size_t data_dim );

    std::vector<double> edgeWeightsLaplace3d( const topology::CombinatorialMap& map,
                                              const VertexPositionsFunc& vertex_position,
                                              const std::vector<Normal>& normals,
                                              const LaplaceEdgeWeights& edge_weights );

} // namespace reparam