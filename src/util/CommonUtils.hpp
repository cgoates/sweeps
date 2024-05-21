#pragma once
#include <Eigen/Dense>

namespace util
{
    bool equals( const double& a, const double& b, const double& tol );

    bool equals( const Eigen::Ref<const Eigen::VectorXd> a,
                 const Eigen::Ref<const Eigen::VectorXd> b,
                 const double& tol );
} // namespace util