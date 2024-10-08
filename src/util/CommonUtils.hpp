#pragma once
#include <Eigen/Dense>

namespace util
{
    bool equals( const double& a, const double& b, const double& tol );

    bool equals( const Eigen::Ref<const Eigen::VectorXd> a,
                 const Eigen::Ref<const Eigen::VectorXd> b,
                 const double& tol );

    double normalizeAngle( double angle );

    bool angleEquals( const double a, const double b, const double tol );

    std::vector<double> linspace( const double left_val, const double right_val, const size_t n_levels );

    std::vector<double> concatenate( const std::vector<double>& first, const std::vector<double>& second );
} // namespace util