#pragma once
#include <SplineSpaceEvaluator.hpp>

namespace fitting
{
    Eigen::MatrixXd leastSquaresFitting(
        eval::SplineSpaceEvaluator& evaler,
        const size_t n_points,
        const size_t spatial_dim,
        const std::function<void( const std::function<void(const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd&)>& )>& point_iterator );
}