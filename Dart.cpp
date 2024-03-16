#include <Dart.hpp>

std::ostream& operator<<( std::ostream& o, const topology::Dart& d )
{
    o << "Dart(" << d.id() << ")";
    return o;
}