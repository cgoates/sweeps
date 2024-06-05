#pragma once
#include <Eigen/Dense>
#include <map>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class Cell;
    class Vertex;
    class Edge;
    class CombinatorialMap;
}

namespace reparam
{
    class Trace;
    using TraceLevelSetIntersection = std::pair<Eigen::Vector3d, size_t>;

    /// @brief Find locations where a trace intersects level sets.
    /// @param trace            The trace that is to be processed.
    /// @param level_set_values The values of the scalar field at which to find the intersections. Should be ordered low to high.
    /// @return A vector of intersection locations.
    std::vector<TraceLevelSetIntersection> levelSetIntersections( const Trace& trace,
                                                                  const topology::CombinatorialMap& traced_cmap,
                                                                  const std::vector<double>& level_set_values );

    /// @brief Find the theta values for a circular Tutte embedding on a level set boundary.
    /// @param level_set    The level set to find the boundaries to map to.
    /// @param level_set_positions  The vertex positions of the level set in physical space.
    /// @param underlying_face_id_of_edge Given an edge on the boundary of level_set, returns the face id
    ///                                   in the underlying topology that level_set is a level set of.
    /// @param intersection The location of the intersection which represents the theta=0 location.
    /// @return A map of the theta values for the boundary vertices.
    std::map<topology::Vertex, double>
        thetaValues( const topology::CombinatorialMap& level_set,
                     const VertexPositionsFunc& level_set_positions,
                     const std::function<size_t( const topology::Edge& )>& underlying_face_id_of_edge,
                     const TraceLevelSetIntersection& intersection );
} // namespace reparam