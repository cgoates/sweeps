#pragma once
#include <Eigen/Core>
#include <set>
#include <VertexPositionsFunc.hpp>
#include <optional>
#include <map>

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
        Uniform,
        EdgeLength,
        AverageAdjacentArea,
        InverseAverageAdjacentArea,
        MeanValue,
        RegularizedInverseLength,
        AverageInverseAdjacentArea
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

    std::pair<Eigen::MatrixX2d, std::vector<double>>
        cutTutteEmbedding( const topology::CutCombinatorialMap& map,
                           const VertexPositionsFunc& vert_positions,
                           const std::vector<std::pair<topology::Vertex, topology::Vertex>>& cut_extremities,
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

    /// @brief Find the boundary constraints for a Tutte embedding on the fundamental domain of a cut cmap.
    /// @param cut_cmap           The cut cmap to find the boundary constraints for.
    /// @param cut_cmap_positions The physical domain positions of the vertices of the cut cmap.
    /// @param n_cuts             The number of cuts made to the original cmap to cut it to a disk.
    /// @param is_cut_extremity   An indicator function for vertices which are the ends of cuts.
    /// @return The boundary vertex constraints to tutte embed cut_cmap to its fundamental domain. Keys are vert ids.
    std::map<size_t, Eigen::Vector2d>
        boundaryConstraints( const topology::CombinatorialMap& cut_cmap,
                             const VertexPositionsFunc& cut_cmap_positions,
                             const size_t n_cuts,
                             const std::function<bool( const topology::Vertex& )>& is_cut_extremity,
                             const std::function<bool( const topology::Vertex& )>& is_start_vert );

} // namespace reparam