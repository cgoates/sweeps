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
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            iterateDartsOfCell( map, c, [&]( const Dart& d ) {
                mDartMarker.mark( d );
                return true;
            } );
        }

        void unmark( const CombinatorialMap& map, const Cell& c )
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            iterateDartsOfCell( map, c, [&]( const Dart& d ) {
                mDartMarker.unmark( d );
                return true;
            } );
        }

        bool isMarked( const Cell& c ) const
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            return mDartMarker.isMarked( c.dart() );
        }

        private:
        GlobalDartMarker mDartMarker;
        uint mCellDim;
    };

    class LocalCellMarker
    {
        public:
        LocalCellMarker( const uint cell_dim ) : mCellDim( cell_dim )
        {}

        void mark( const CombinatorialMap& map, const Cell& c )
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            iterateDartsOfCell( map, c, [&]( const Dart& d ) {
                mDartMarker.mark( d );
                return true;
            } );
        }

        bool isMarked( const Cell& c ) const
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            return mDartMarker.isMarked( c.dart() );
        }

        private:
        LocalDartMarker mDartMarker;
        uint mCellDim;
    };

    class IndexingCellMarker
    {
        public:
        IndexingCellMarker( const IndexingFunc& indexing, const uint cell_dim ) : mIndexing( indexing ), mCellDim( cell_dim )
        {}

        void mark( const CombinatorialMap&, const Cell& c )
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            mMarkedCellIds.push_back( mIndexing( c ) );
        }

        bool isMarked( const Cell& c ) const
        {
            if( c.dim() != mCellDim ) throw std::runtime_error( "Bad cell dimension!" );
            return std::find( mMarkedCellIds.begin(), mMarkedCellIds.end(), mIndexing( c ) ) != mMarkedCellIds.end();
        }

        private:
        IndexingFunc mIndexing;
        GrowableVector<size_t, 300> mMarkedCellIds;
        uint mCellDim;
    };
} // namespace topology