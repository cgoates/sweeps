#pragma once
#include <Dart.hpp>

using uint = unsigned int;

namespace topology
{
    class Cell
    {
        public:
        Cell( const Dart& d, const uint dim ) : mDart( d ), mDim( dim ) {}
        Cell() {}

        const Dart& dart() const { return mDart; }
        uint dim() const { return mDim; }

        bool operator<( const Cell& o ) const { return dim() < o.dim() or dart() < o.dart(); }
        bool operator==( const Cell& o ) const { return dim() == o.dim() and dart() == o.dart(); }

        private:
        Dart mDart;
        uint mDim;
    };

    class Vertex : public Cell
    {
        public:
        Vertex( const Dart& d ) : Cell( d, 0 ) {}
        Vertex( const Cell& c );
        Vertex() {}
    };

    class Edge : public Cell
    {
        public:
        Edge( const Dart& d ) : Cell( d, 1 ) {}
        Edge( const Cell& c );
        Edge() {}
    };

    class Face : public Cell
    {
        public:
        Face( const Dart& d ) : Cell( d, 2 ) {}
        Face( const Cell& c );
        Face() {}
    };

    class Volume : public Cell
    {
        public:
        Volume( const Dart& d ) : Cell( d, 3 ) {}
        Volume( const Cell& c );
        Volume() {}
    };

    std::ostream& operator<<( std::ostream& o, const topology::Cell& d );
} // namespace topology
