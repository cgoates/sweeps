#pragma once
#include <Eigen/Core>
#include <numeric>

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

    template<typename T>
    std::vector<T> concatenate( const std::initializer_list<const std::vector<T>>& vecs )
    {
        std::vector<T> out;
        const size_t n = std::transform_reduce(
            vecs.begin(), vecs.end(), 0, std::plus<>(), []( const auto& v ) { return v.size(); } );
        out.reserve( n );

        for( const std::vector<T>& v : vecs )
            out.insert( out.end(), v.begin(), v.end() );

        return out;
    }

    std::vector<Eigen::Vector2d> regularNGonVertices( const size_t n_sides );

    /// \brief Generates n_points points in a regular n_sides-gon.
    /// These points are randomly distributed, but are consistent between calls.
    std::vector<Eigen::Vector2d> generatePointsInPolygon( const size_t n_points, const size_t n_sides );

    /// \brief Generates n_points randomly distributed in a circle of radius 1 centered at the origin.
    ///
    std::vector<Eigen::Vector2d> generatePointsInCircle( const size_t n_points );
} // namespace util