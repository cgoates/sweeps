#pragma once
#include <VertexPositionsFunc.hpp>

namespace basis
{
    class SplineSpace;
}

namespace fitting
{
    /// @brief Constructs control points for a linear spline space from vertex positions.
    /// @param ss The spline space for which to construct control points. Must be of order 1 and be locally linearly independent.
    /// @param positions A function that returns the position of a vertex.
    /// @return A matrix of control points, where each row corresponds to a function in the spline space.
    Eigen::MatrixXd linearControlPointsFromVertexPositions( const basis::SplineSpace& ss, const VertexPositionsFunc& positions );
}