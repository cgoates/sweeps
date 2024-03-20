#include <Cell.hpp>


std::ostream& operator<<( std::ostream& o, const topology::Cell& c )
{
    switch( c.dim() )
    {
        case 0:
            o << "Vertex( ";
            break;
        case 1:
            o << "Edge( ";
            break;
        case 2:
            o << "Face( ";
            break;
        case 3:
            o << "Volume( ";
            break;
        default:
            o << c.dim() << "-cell( ";
    }
    o << c.dart() << " )";
    return o;
}