#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalDartMarker.hpp>
#include <GlobalCellMarker.hpp>
#include <queue>
#include <Logging.hpp>

using namespace topology;

CombinatorialMapBoundary::CombinatorialMapBoundary( const CombinatorialMap& map ) : 
    CombinatorialMapBoundary( map, boundaryComponentDarts( map ) )
{}

CombinatorialMapBoundary::CombinatorialMapBoundary( const CombinatorialMap& map, const std::set<Dart>& boundary_component_darts ) :
    mInteriorMap( map ), mStartDarts( boundary_component_darts )
{
    if( map.dim() < 2 or map.dim() > 3 )
        throw( "Unsupported map dimension for creating a boundary" );
}

std::optional<Dart> CombinatorialMapBoundary::phi( const int i, const Dart& d ) const
{
    const int effective_i = std::abs( i ) == 1 ? -i : i;
    if( i < (int)dim() ) return topology::phi( mInteriorMap, effective_i, d );
    else if( i > (int)dim() ) return std::nullopt;
    else
    {
        const std::vector<int> phi_ops = { effective_i, std::abs( i ) + 1 };
        Dart curr_dart = d;
        std::optional<Dart> next_dart = d;
        while( next_dart.has_value() )
        {
            curr_dart = next_dart.value();
            next_dart = topology::phi( mInteriorMap, phi_ops, curr_dart );
        }
        return topology::phi( mInteriorMap, effective_i, curr_dart );
    }
}

Dart::IndexType CombinatorialMapBoundary::maxDartId() const
{
    return mInteriorMap.maxDartId();
}

uint CombinatorialMapBoundary::dim() const
{
    return mInteriorMap.dim() - 1;
}

bool CombinatorialMapBoundary::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    std::queue<Dart> to_flood;
    GlobalDartMarker m( *this );

    for( const Dart& d : mStartDarts )
    {
        to_flood.push( d );
        m.mark( d );
    }

    while( not to_flood.empty() )
    {
        const Dart d = to_flood.front();
        to_flood.pop();
        if( not callback( d ) ) return false;
        for( int i : { 1, 2 } )
        {
            const auto a = topology::phi( *this, i, d );
            if( a.has_value() and not m.isMarked( a.value() ) )
            {
                to_flood.push( a.value() );
                m.mark( a.value() );
            }
        }
    }

    return true;
}

bool CombinatorialMapBoundary::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
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