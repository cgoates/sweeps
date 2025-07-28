#pragma once
#include <SplineSpaceEvaluator.hpp>

namespace basis
{
    class MultiPatchSplineSpace;
}

namespace fitting
{
    Eigen::MatrixXd leastSquaresFitting(
        eval::SplineSpaceEvaluator& evaler,
        const size_t n_points,
        const size_t spatial_dim,
        const std::function<void( const std::function<void(const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd&)>& )>& point_iterator );

    /// @brief Fits a spline space to a manifold defined over the same parametric atlas.
    /// @param manifold_ss   The spline space defining the manifold.
    /// @param manifold_cpts The control points defining the manifold.
    /// @param new_ss        The spline space to fit to the manifold.
    /// @return The control points of the new spline space fitted to the manifold.
    Eigen::MatrixXd fitToManifold( const basis::MultiPatchSplineSpace& manifold_ss,
                                   const Eigen::MatrixXd& manifold_cpts,
                                   const basis::MultiPatchSplineSpace& new_ss );
}