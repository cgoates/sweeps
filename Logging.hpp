#pragma once

#define LOG( COND ) if( COND ) std::cout

template<typename T>
std::ostream& operator<<( std::ostream& o, const std::vector<T>& v )
{
    if( v.size() == 0 ) o << "{}";
    else
    {    
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << *it << ", ";
        o << v.back() << " }";
    }
    return o;
}

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

template<typename T>
std::ostream& operator<<( std::ostream& o, const std::set<T>& v )
{
    if( v.size() == 0 ) o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end(); it++ ) o << *it << ", ";
        o << " }";
    }
    return o;
}