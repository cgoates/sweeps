#pragma once
#include <Eigen/Dense>
#include <map>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class Cell;
    class Vertex;
    class LevelSetCMap;
}

namespace reparam
{
    class Trace;
    using TraceLevelSetIntersection = std::pair<Eigen::Vector3d, topology::Cell>;

    /// @brief Find locations where a trace intersects level sets.
    /// @param trace            The trace that is to be processed.
    /// @param level_set_values The values of the scalar field at which to find the intersections. Should be ordered low to high.
    /// @return A vector of intersection locations.
    std::vector<TraceLevelSetIntersection> levelSetIntersections( const Trace& trace,
                                                                  const std::vector<double>& level_set_values );

    /// @brief Find the theta values for a circular Tutte embedding on a level set boundary.
    /// @param level_set    The level set to find the boundaries to map to.
    /// @param intersection The location of the intersection which represents the theta=0 location.
    /// @return A map of the theta values for the boundary vertices.
    /// FIXME: This should also be able to take a source or target surface
    std::map<topology::Vertex, double> thetaValues( const topology::LevelSetCMap& level_set,
                                                    const VertexPositionsFunc& level_set_positions,
                                                    const TraceLevelSetIntersection& intersection );
} // namespace reparam