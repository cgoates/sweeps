#pragma once
#include <Eigen/Core>
#include <set>
#include <VertexPositionsFunc.hpp>
#include <optional>

class VertexId;
class Normal;

namespace topology
{
    class TetMeshCombinatorialMap;
    class CutCombinatorialMap;
    class CombinatorialMap;
    class Edge;
}

namespace reparam
{
    enum class Laplace3dEdgeWeights
    {
        Cotangent,
        InverseLength,
        VoronoiDual,
        BarycentricDual,
        Uniform
    };

    enum class Laplace2dEdgeWeights
    {
        Cotangent,
        InverseLength,
        BarycentricDual,
        Uniform
    };

    Eigen::VectorXd sweepEmbedding( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals,
                                    const Laplace3dEdgeWeights& edge_weights = Laplace3dEdgeWeights::BarycentricDual );

    Eigen::MatrixX2d
        tutteEmbedding( const topology::CombinatorialMap& map,
                        const VertexPositionsFunc& vert_positions,
                        const std::function<std::optional<Eigen::Vector2d>( const topology::Vertex& )>& constraints,
                        const Laplace2dEdgeWeights& edge_weights_type = Laplace2dEdgeWeights::InverseLength );

    Eigen::MatrixX2d tutteOrbifoldEmbedding( const topology::CutCombinatorialMap& map,
                                             const VertexPositionsFunc& vert_positions,
                                             const std::array<topology::Vertex, 3>& cone_vertices,
                                             const Laplace2dEdgeWeights& edge_weights_type = Laplace2dEdgeWeights::InverseLength );

    Eigen::MatrixXd
        solveLaplaceSparse( const topology::CombinatorialMap& map,
                            const std::function<double( const topology::Edge& )>& edge_weights,
                            const std::function<std::optional<Eigen::VectorXd>( const topology::Vertex& )>& constraints,
                            const size_t n_constrained_verts,
                            const size_t data_dim );

    std::vector<double> edgeWeightsLaplace3d( const topology::CombinatorialMap& map,
                                              const VertexPositionsFunc& vertex_position,
                                              const std::vector<Normal>& normals,
                                              const Laplace3dEdgeWeights& edge_weights );

} // namespace reparam