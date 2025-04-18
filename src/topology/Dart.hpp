#pragma once
#include <cstdint>
#include <iosfwd>

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

        bool operator!=( const Dart& o ) const { return id() != o.id(); }

        private:
        IndexType mIndex;
    };

    std::ostream& operator<<( std::ostream& o, const topology::Dart& d );
}; // namespace topology