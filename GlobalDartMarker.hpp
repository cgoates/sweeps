#pragma once
#include <CombinatorialMap.hpp>
#include <vector>
#include <set>

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

        void mark( const Dart& d ) { mMarkedDarts.insert( d ); }

        bool isMarked( const Dart& d ) const { return std::find( mMarkedDarts.begin(), mMarkedDarts.end(), d ) != mMarkedDarts.end(); }

        private:
        std::set<Dart> mMarkedDarts;
    };
} // namespace topology