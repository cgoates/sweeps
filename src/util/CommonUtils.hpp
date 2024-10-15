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

    template<typename T>
    std::vector<T> concatenate( const std::vector<T>& first, const std::vector<T>& second )
    {
        std::vector<T> out;
        out.reserve( first.size() + second.size() );

        out.insert( out.end(), first.begin(), first.end() );
        out.insert( out.end(), second.begin(), second.end() );

        return out;
    }
} // namespace util