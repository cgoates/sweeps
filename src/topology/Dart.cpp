#include <Dart.hpp>
#include <iostream>

namespace topology{
    std::ostream& operator<<( std::ostream& o, const topology::Dart& d )
    {
        o << "Dart(" << d.id() << ")";
        return o;
    }
}