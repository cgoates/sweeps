#include <CustomCombinatorialMap.hpp>
#include <GlobalCellMarker.hpp>

namespace topology
{
    CustomCombinatorialMap::CustomCombinatorialMap( const size_t n_darts,
                            const uint dim,
                            const SmallVector<std::vector<std::optional<Dart::IndexType>>, 3>& phis,
                            const std::map<size_t, std::vector<size_t>>& cell_ids ) :
        mCellIds( cell_ids ),
        mDim( dim )
    {
        mPhi1s = phis.at( 0 );
        if( phis.size() > 1 )
        {
            mPhi2s = phis.at( 1 );
            if( phis.size() > 2 )
                mPhi3s = phis.at( 2 );
        }

        mPhi_1s = std::vector<std::optional<Dart::IndexType>>( n_darts, std::nullopt );
        for( size_t i = 0; i < n_darts; i++ )
        {
            if( mPhi1s.at( i ).has_value() )
            {
                mPhi_1s.at( mPhi1s.at( i ).value() ) = i;
            }
        }
    }

    std::optional<Dart> CustomCombinatorialMap::phi( const int i, const Dart& d ) const
    {
        const std::vector<std::optional<Dart::IndexType>>& phi = ( i == -1 ) ? mPhi_1s : ( i == 1 ) ? mPhi1s : ( i == 2 ) ? mPhi2s : mPhi3s;
        return phi.at( d.id() ).transform( [&]( const Dart::IndexType& id ) { return Dart( id ); } );
    }
    Dart::IndexType CustomCombinatorialMap::maxDartId() const
    {
        return mPhi1s.size() - 1;
    }
    uint CustomCombinatorialMap::dim() const { return mDim; }
    bool CustomCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
    {
        for( size_t i = 0; i < mPhi1s.size(); i++ )
        {
            if( not callback( Dart( i ) ) ) return false;
        }
        return true;
    }
    bool CustomCombinatorialMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
    {
        GlobalCellMarker m( *this, cell_dim );
        return iterateDartsWhile( [&]( const Dart& d ) {
            const Cell c( d, cell_dim );
            if( m.isMarked( c ) ) return true;
            if( not callback( c ) ) return false;
            m.mark( *this, c );
            return true;
        } );
    }
    std::optional<IndexingFunc> CustomCombinatorialMap::indexing( const uint cell_dim ) const
    {
        const auto it = mCellIds.find( cell_dim );
        if( it == mCellIds.end() ) return std::nullopt;
        return [&ids = it->second]( const Cell& c ) -> size_t {
            return ids.at( c.dart().id() );
        };
    }

    std::optional<size_t> CustomCombinatorialMap::cellCount( const uint ) const
    {
        return std::nullopt;
    }
} // namespace topology