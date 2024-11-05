#include <CommonUtils.hpp>
#include <numbers>

namespace util
{

    bool equals( const double& a, const double& b, const double& tol )
    {
        return std::abs( a - b ) < tol;
    }

    bool equals( const Eigen::Ref<const Eigen::VectorXd> a,
                 const Eigen::Ref<const Eigen::VectorXd> b,
                 const double& tol )
    {
        if( a.rows() != b.rows() ) return false;
        for( Eigen::Index i = 0; i < b.rows(); i++ )
        {
            if( not equals( a( i ), b( i ), tol ) ) return false;
        }
        return true;
    }

    double normalizeAngle( double angle )
    {
        while( angle < -1 * std::numbers::pi ) angle += 2 * std::numbers::pi;
        while( angle > std::numbers::pi ) angle -= 2 * std::numbers::pi;
        return angle;
    }

    bool angleEquals( const double a, const double b, const double tol )
    {
        return equals( normalizeAngle( a ), normalizeAngle( b ), tol );
    }

    std::vector<double> linspace( const double left_val, const double right_val, const size_t n_levels )
    {
        std::vector<double> out;
        out.reserve( n_levels );
        const double diff = ( right_val - left_val ) / ( n_levels - 1 );
        for( size_t i = 0; i < n_levels; i++ )
        {
            out.push_back( left_val + i * diff );
        }
        return out;
    }

    std::vector<Eigen::Vector2d> regularNGonVertices( const size_t n_sides )
    {
        std::vector<Eigen::Vector2d> out;
        out.reserve( n_sides + 1 ); // Add endpoint twice to avoid having to worry about wraparound.
        for( size_t i = 0; i <= n_sides; i++ )
        {
            const double theta = 2 * i * std::numbers::pi / n_sides;
            out.push_back( { cos( theta ), sin( theta ) } );
        }
        return out;
    }
} // namespace util