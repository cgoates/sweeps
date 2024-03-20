#include <Logging.hpp>

std::ostream& operator<<( std::ostream& o, const Eigen::Triplet<double>& t )
{
    o << "{" << t.row() << ", " << t.col() << ": " << t.value() << "}";
    return o;
}

std::ostream& operator<<( std::ostream& o, const std::vector<Eigen::Vector3d>& v )
{
    if( v.size() == 0 )
        o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << it->transpose() << ",\n  ";
        o << v.back().transpose() << " }";
    }
    return o;
}

bool equals( const double& a, const double& b, const double& tol )
{
    return std::abs( a - b ) < tol;
}

bool equals( const Eigen::Ref<const Eigen::VectorXd> a, const Eigen::Ref<const Eigen::VectorXd> b, const double& tol )
{
    if( a.rows() != b.rows() ) return false;
    for( Eigen::Index i = 0; i < b.rows(); i++ )
    {
        if( not equals( a( i ), b( i ), tol ) ) return false;
    }
    return true;
}