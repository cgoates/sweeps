#pragma once
#include <CombinatorialMap.hpp>
#include <vector>
#include <SmallVector.hpp>

namespace topology
{
    class GlobalDartMarker
    {
        public:
        GlobalDartMarker( const CombinatorialMap& map ) : mMarkedDarts( map.maxDartId() + 1, false ) {}

        void mark( const Dart& d ) { mMarkedDarts.at( d.id() ) = true; }

        bool isMarked( const Dart& d ) const { return mMarkedDarts.at( d.id() ); }

        private:
        std::vector<bool> mMarkedDarts;
    };

    class LocalDartMarker
    {
        public:
        LocalDartMarker() {}

        void mark( const Dart& d ) { mMarkedDarts.push_back( d ); }

        bool isMarked( const Dart& d ) const { return std::find( mMarkedDarts.begin(), mMarkedDarts.end(), d ) != mMarkedDarts.end(); }

        private:
        // FIXME: I need a growable vector that uses a small stack allocation, but allocates on the heap if needed.
        //  I keep having to up the number on this, but that removes some of the benefit for the majority of the cases
        //  where it doesn't need as large of a vector.
        SmallVector<Dart, 1200> mMarkedDarts;
    };
} // namespace topology