#include<Logging.hpp>

std::ostream& operator<<( std::ostream& o, const std::vector<Eigen::Vector3d>& v )
{
    if( v.size() == 0 ) o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << it->transpose() << ",\n  ";
        o << v.back().transpose() << " }";
    }
    return o;
}