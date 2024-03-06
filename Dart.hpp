#pragma once
#include <cstdint>

namespace topology
{
    class Dart
    {
        public:
        using IndexType = uint64_t;

        Dart( const IndexType& id ) : mIndex( id ) {}

        IndexType id() const { return mIndex; }

        bool operator<( const Dart& o ) const { return id() < o.id(); }

        private:
        IndexType mIndex;
    };
}; // namespace topology