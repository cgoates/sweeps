#pragma once
#include <CombinatorialMap.hpp>
#include <GlobalDartMarker.hpp>
#include <CombinatorialMapMethods.hpp>

namespace topology
{
    class GlobalCellMarker
    {
        public:
        GlobalCellMarker( const CombinatorialMap& map, const uint cell_dim ) : mDartMarker( map ), mCellDim( cell_dim )
        {}

        void mark( const CombinatorialMap& map, const Cell& c )
        {
            if( c.dim() != mCellDim ) throw( "Bad cell dimension!" );
            iterateDartsOfCell( map, c, [&]( const Dart& d ) {
                mDartMarker.mark( d );
                return true;
            } );
        }

        bool isMarked( const Cell& c ) const
        {
            if( c.dim() != mCellDim ) throw( "Bad cell dimension!" );
            return mDartMarker.isMarked( c.dart() );
        }

        private:
        GlobalDartMarker mDartMarker;
        uint mCellDim;
    };
} // namespace topology