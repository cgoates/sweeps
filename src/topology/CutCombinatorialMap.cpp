#include <CutCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>

using namespace topology;

CutCombinatorialMap::CutCombinatorialMap( const CombinatorialMap& base, const std::set<topology::Cell>& interfaces_to_disconnect )
    : mBaseMap( base ), mNumCuts( interfaces_to_disconnect.size() )
{
    for( const auto& c : interfaces_to_disconnect )
    {
        if( c.dim() != dim() - 1 )
            throw std::invalid_argument( "Can only cut a cmap along (dim-1)-cells" );

        iterateDartsOfCell( base, c, [&]( const Dart& d ){
            mNoPhiDimDarts.insert( d );
            return true;
        } );
    }

    // Build vertex indexing
    const auto maybe_ids = mBaseMap.indexing( 0 );
    if( maybe_ids.has_value() )
    {
        const auto& vert_ids = maybe_ids.value();
        GlobalCellMarker m( *this, 0 );
        size_t next_id = 0;
        topology::iterateCellsWhile( mBaseMap, 0, [&]( const Vertex& v ){
            m.mark( *this, v );
            next_id = std::max( next_id, vert_ids( v ) );
            return true;
        } );
        next_id++;
        for( const Dart& d : mNoPhiDimDarts )
        {
            iterateAdjacentCells( *this, Cell( d, dim() - 1 ), 0, [&]( const Vertex& v ) {  
                if( not m.isMarked( v ) )
                {
                    m.mark( *this, v );
                    iterateDartsOfCell( *this, v, [&]( const Vertex& equivalent_v ){
                        mAdditionalVertexIds.emplace( equivalent_v, next_id );
                        return true;
                    } );
                    next_id++;
                }
                return true;
            } );
        }
    }
}

std::optional<Dart> CutCombinatorialMap::phi( const int i, const Dart& d ) const
{
    if( i == (int)dim() and mNoPhiDimDarts.contains( d ) ) return std::nullopt;
    else
        return topology::phi( mBaseMap, i, d );
}

bool CutCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    return topology::iterateDartsWhile( mBaseMap, callback );
}

bool CutCombinatorialMap::iterateCellsWhile( const uint cell_dim,
                                               const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == dim() )
    {
        return topology::iterateCellsWhile( mBaseMap, cell_dim, callback );
    }
    else
    {
        GlobalCellMarker m( *this, cell_dim );
        const bool keep_iterating = topology::iterateCellsWhile( mBaseMap, cell_dim, [&]( const Cell& c ){
            m.mark( *this, c );
            return callback( c );
        } );
        if( keep_iterating )
        {
            for( const Dart& d : mNoPhiDimDarts )
            {
                const bool keep_iterating2 = iterateAdjacentCells( *this, Cell( d, dim() - 1 ), cell_dim, [&]( const Cell& c ) {
                    if( not m.isMarked( c ) )
                    {
                        m.mark( *this, c );
                        if( not callback( c ) ) return false;
                    }
                    return true;
                } );
                if( not keep_iterating2 ) return false;
            }
            return true;
        }
        else
            return false;
    }
}

std::optional<IndexingFunc> CutCombinatorialMap::indexing( const uint cell_dim ) const
{
    if( cell_dim == dim() )
        return mBaseMap.indexing( cell_dim );
    else if( cell_dim == 0 )
    {
        return mBaseMap.indexing( cell_dim ).transform( [&]( const auto& base_ids ) -> IndexingFunc {
            return [this,base_ids]( const Vertex& v ){
                const auto it = mAdditionalVertexIds.find( v );
                if( it == mAdditionalVertexIds.end() )
                    return base_ids( v );
                else
                    return it->second;
            };
        } );
    }

    return std::nullopt;
}

std::optional<size_t> CutCombinatorialMap::cellCount( const uint cell_dim ) const
{
    if( cell_dim == dim() ) return topology::cellCount( mBaseMap, cell_dim );
    else if( cell_dim == dim() - 1 ) return topology::cellCount( mBaseMap, cell_dim ) + mNumCuts;
    else return std::nullopt;
}