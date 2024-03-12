#pragma once
#include <cstdint>
#include <string> // For hash

namespace topology
{
    class Dart
    {
        public:
        using IndexType = uint64_t;

        Dart( const IndexType& id ) : mIndex( id ) {}
        Dart() {}

        IndexType id() const { return mIndex; }

        bool operator<( const Dart& o ) const { return id() < o.id(); }

        bool operator==( const Dart& o ) const { return id() == o.id(); }

        private:
        IndexType mIndex;
    };
}; // namespace topology

template <>
struct std::hash<topology::Dart>
{
    std::size_t operator()(const topology::Dart& d) const
    {
        return hash<topology::Dart::IndexType>()( d.id() );
    }
};