#pragma once
#include <CombinatorialMap.hpp>
#include <vector>

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
} // namespace topology