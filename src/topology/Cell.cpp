#include <Cell.hpp>
#include <iostream>

namespace topology
{
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

Vertex::Vertex( const Cell& c ) : Cell( c )
{
    if( c.dim() != 0 )
        throw std::runtime_error( "Cannot convert cell with dimension " + std::to_string( c.dim() ) + " to Vertex" );
}

Edge::Edge( const Cell& c ) : Cell( c )
{
    if( c.dim() != 1 )
        throw std::runtime_error( "Cannot convert cell with dimension " + std::to_string( c.dim() ) + " to Edge" );
}

Face::Face( const Cell& c ) : Cell( c )
{
    if( c.dim() != 2 )
        throw std::runtime_error( "Cannot convert cell with dimension " + std::to_string( c.dim() ) + " to Face" );
}

Volume::Volume( const Cell& c ) : Cell( c )
{
    if( c.dim() != 3 )
        throw std::runtime_error( "Cannot convert cell with dimension " + std::to_string( c.dim() ) + " to Volume" );
}

}