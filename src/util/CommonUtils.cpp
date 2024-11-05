#include <CommonUtils.hpp>
#include <numbers>
#include <random>

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

    std::vector<Eigen::Vector2d> generatePointsInPolygon( const size_t n_points, const size_t n_sides )
    {

        const auto is_inside_polygon = []( const Eigen::Vector2d& point, const std::vector<Eigen::Vector2d>& vertices ) {
            for( size_t i = 0; i < vertices.size() - 1; ++i )
            {
                Eigen::Vector2d edge = vertices.at( i + 1 ) - vertices.at( i );
                Eigen::Vector2d to_point = point - vertices.at( i );
                if( edge( 0 ) * to_point( 1 ) - edge( 1 ) * to_point( 0 ) < 0 )
                {
                    return false;
                }
            }
            return true;
        };
        std::vector<Eigen::Vector2d> points;
        points.reserve( n_points );

        const auto vertices = util::regularNGonVertices( n_sides );

        // Use a seed to get the same points every time
        const size_t seed = 0;
        std::mt19937 gen( seed );
        std::uniform_real_distribution<> dis( -1.0, 1.0 );

        while( points.size() < n_points )
        {
            Eigen::Vector2d point( dis( gen ), dis( gen ) );
            if( is_inside_polygon( point, vertices ) )
            {
                points.push_back( point );
            }
        }

        return points;
    }

    std::vector<Eigen::Vector2d> generatePointsInCircle( const size_t n_points )
    {

        const auto is_inside_polygon = []( const Eigen::Vector2d& point ) { return point.norm() < 1; };
        std::vector<Eigen::Vector2d> points;
        points.reserve( n_points );

        // Use a seed to get the same points every time
        const size_t seed = 0;
        std::mt19937 gen( seed );
        std::uniform_real_distribution<> dis( -1.0, 1.0 );

        while( points.size() < n_points )
        {
            Eigen::Vector2d point( dis( gen ), dis( gen ) );
            if( is_inside_polygon( point ) )
            {
                points.push_back( point );
            }
        }

        return points;
    }
} // namespace util