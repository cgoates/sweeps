#pragma once
#include <array>
#include <span>
#include <iostream>
#include <ecpp/static_vector.hpp>

template<typename T, size_t MAX_SIZE>
using SmallVector = ecpp::static_vector<T, MAX_SIZE>;

namespace ecpp
{
template <typename T, size_t N> std::ostream& operator<<( std::ostream& o, const ecpp::static_vector<T, N>& v )
{
    if( v.size() == 0 )
        o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << *it << ", ";
        o << *( v.end() - 1 ) << " }";
    }
    return o;
}
}