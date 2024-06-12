#include <CombinatorialMapRestriction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Simplex.hpp>

using namespace topology;

CombinatorialMapRestriction::CombinatorialMapRestriction( const CombinatorialMap& map,
                                                          const std::function<bool( const Cell& )>& restriction_func,
                                                          const bool reindex_verts )
    : mUnrestrictedMap( map ), mIncludedDarts( map )
{
    topology::iterateCellsWhile( map, map.dim(), [&]( const Cell& elem ) {
        if( restriction_func( elem ) )
        {
            iterateDartsOfCell( map, elem, [&]( const Dart& d ) {
                mIncludedDarts.mark( d );
                return true;
            } );
        }
        return true;
    } );

    if( reindex_verts )
    {
        const auto underlying_ids = map.indexing( 0 );
        if( underlying_ids.has_value() )
        {
            mVertexIds.emplace( std::map<size_t, size_t>() );
            size_t vert_ii = 0;
            iterateCellsWhile( 0, [&]( const topology::Vertex& v ) {
                mVertexIds->insert( { underlying_ids.value()( v ), vert_ii++ } );
                return true;
            } );
        }

    }
}

std::optional<Dart> CombinatorialMapRestriction::phi( const int i, const Dart& d ) const
{
    const auto candidate = topology::phi( mUnrestrictedMap, i, d );
    if( std::abs( i ) != (int)dim() or not candidate.has_value() or mIncludedDarts.isMarked( candidate.value() ) )
    {
        return candidate;
    }
    return std::nullopt;
}

bool CombinatorialMapRestriction::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    return topology::iterateDartsWhile( mUnrestrictedMap, [&]( const Dart& d ) {
        if( mIncludedDarts.isMarked( d ) )
        {
            if( not callback( d ) ) return false;
        }
        return true;
    } );
}

bool CombinatorialMapRestriction::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    return topology::iterateCellsWhile( mUnrestrictedMap, cell_dim, [&]( const Cell& c ) {
        if( mIncludedDarts.isMarked( c.dart() ) )
        {
            if( not callback( c ) ) return false;
        }
        else
        {
            bool continue_iterating = true;
            iterateDartsOfCell( mUnrestrictedMap, c, [&]( const Dart& d ) {
                if( mIncludedDarts.isMarked( d ) )
                {
                    continue_iterating = callback( Cell( d, cell_dim ) );
                    return false;
                }
                return true;
            } );
            if( not continue_iterating ) return false;
        }
        return true;
    } );
}

std::optional<IndexingFunc> CombinatorialMapRestriction::indexing( const uint cell_dim ) const
{
    if( cell_dim == 0 and mVertexIds )
    {
        return mUnrestrictedMap.indexing( cell_dim ).transform( [&]( const auto& underlying_ids ) -> IndexingFunc {
            return [underlying_ids,this]( const topology::Vertex& v ){
                return mVertexIds->at( underlying_ids( v ) );
            };
        } );
    }
    return mUnrestrictedMap.indexing( cell_dim );
}