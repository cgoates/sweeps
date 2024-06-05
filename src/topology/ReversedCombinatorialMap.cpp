#include <ReversedCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalDartMarker.hpp>
#include <GlobalCellMarker.hpp>
#include <queue>
#include <Logging.hpp>
#include <Simplex.hpp>

using namespace topology;

ReversedCombinatorialMap::ReversedCombinatorialMap( const CombinatorialMap& map ) : 
    mUnderlyingMap( map )
{
    if( map.dim() > 2 )
        throw std::runtime_error( "Unsupported map dimension for reversing" );
}

std::optional<Dart> ReversedCombinatorialMap::phi( const int i, const Dart& d ) const
{
    const int effective_i = std::abs( i ) == 1 ? -i : i;
    return topology::phi( mUnderlyingMap, effective_i, d );
}

Dart::IndexType ReversedCombinatorialMap::maxDartId() const
{
    return mUnderlyingMap.maxDartId();
}

uint ReversedCombinatorialMap::dim() const
{
    return mUnderlyingMap.dim();
}

bool ReversedCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    return topology::iterateDartsWhile( mUnderlyingMap, callback );
}

bool ReversedCombinatorialMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 0 )
    {
        GlobalCellMarker m( *this, cell_dim );
        return iterateDartsWhile( [&]( const Dart& d ){
            const Cell c( d, cell_dim );
            if( not m.isMarked( c ) )
            {
                m.mark( *this, c );
                if( not callback( c ) ) return false;
            }
            return true;
        } );
    }
    else
    {
        return topology::iterateCellsWhile( mUnderlyingMap, cell_dim, callback );
    }
}

topology::Cell ReversedCombinatorialMap::toUnderlyingCell( const topology::Cell& c ) const
{
    if( c.dim() > 0 ) return c;
    else return Vertex( topology::phi( mUnderlyingMap, 1, c.dart() ).value() );
}

std::optional<IndexingFunc> ReversedCombinatorialMap::indexing( const uint cell_dim ) const
{
    return mUnderlyingMap.indexing( cell_dim ).and_then( [&]( const IndexingFunc underlying_func ) -> std::optional<IndexingFunc> {
        return [underlying_func,this]( const Cell& c ) {
            return underlying_func( toUnderlyingCell( c ) );
        };
    } );
}

namespace topology
{
    VertexPositionsFunc reversedVertexPositions( const ReversedCombinatorialMap& bdry,
                                                 const VertexPositionsFunc& underlying_positions )
    {
        return [underlying_positions,&bdry]( const Vertex& v ){
            return underlying_positions( bdry.toUnderlyingCell( v ) );
        };
    }
}