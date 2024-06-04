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
} // namespace util