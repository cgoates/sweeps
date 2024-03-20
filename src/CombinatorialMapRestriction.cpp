#include <CombinatorialMapRestriction.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Simplex.hpp>

using namespace topology;

CombinatorialMapRestriction::CombinatorialMapRestriction( const CombinatorialMap& map,
                                                          const std::function<bool( const Cell& )>& restriction_func )
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
            const bool continue_iterating = iterateDartsOfCell( mUnrestrictedMap, c, [&]( const Dart& d ) {
                if( mIncludedDarts.isMarked( d ) )
                {
                    if( not callback( Cell( d, cell_dim ) ) ) return false;
                }
                return true;
            } );
            if( not continue_iterating ) return false;
        }
        return true;
    } );
}

VertexId CombinatorialMapRestriction::vertexId( const Vertex& v ) const
{
    return mUnrestrictedMap.vertexId( v );
}