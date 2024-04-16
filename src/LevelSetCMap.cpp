#include <LevelSetCMap.hpp>
#include <queue>
#include <Simplex.hpp>

using namespace topology;

Edge findIntersectingEdge( const CombinatorialMap& cmap,
                           const std::function<double( const Vertex& )>& func,
                           const double value )
{
    std::optional<Edge> out = std::nullopt;
    iterateCellsWhile( cmap, 1, [&]( const Edge& e ) {
        const double left_val = func( Vertex( e.dart() ) );
        const double right_val = func( Vertex( phi( cmap, 1, e.dart() ).value() ) );
        if( ( left_val >= value ) != ( right_val >= value ) )
        {
            out.emplace( e );
            return false;
        }
        return true;
    } );
    if( not out.has_value() ) throw std::runtime_error( "No intersecting edge!" );
    return out.value();
}

LevelSetCMap::LevelSetCMap( const CombinatorialMap& base,
                            const std::function<double( const Vertex& )>& func,
                            const double value )
    : LevelSetCMap( base, func, value, findIntersectingEdge( base, func, value ) )
{}

LevelSetCMap::LevelSetCMap( const CombinatorialMap& base,
                            const std::function<double( const Vertex& )>& func,
                            const double value,
                            const Edge& one_intersected_edge )
    : mUnderlyingMap( base ), mIntersectedEdges( base, 1 )
{
    // Fill in mIntersectedEdges and mIntersectionPositions.
    const uint base_dim = base.dim();
    std::queue<Edge> edges_to_process;
    edges_to_process.push( one_intersected_edge );

    while( not edges_to_process.empty() )
    {
        const Edge e = edges_to_process.front();
        edges_to_process.pop();
        iterateAdjacentCells( base, e, base_dim, [&]( const Cell& elem ) {
            iterateAdjacentCells( base, elem, 1, [&]( const Edge& e_adj ) {
                if( not mIntersectedEdges.isMarked( e_adj ) )
                {
                    const double left_val = func( Vertex( e_adj.dart() ) );
                    const double right_val = func( Vertex( topology::phi( base, 1, e_adj.dart() ).value() ) );
                    if( ( left_val >= value ) != ( right_val >= value ) )
                    {
                        edges_to_process.push( e_adj );
                        mIntersectedEdges.mark( base, e_adj );
                        if( left_val >= value )
                        {
                            const double barycentric = ( value - left_val ) / ( right_val - left_val );
                            mIntersectionPositions.emplace( Vertex( e_adj.dart() ), barycentric );
                        }
                        else
                        {
                            const auto maybe_phi2 = topology::phi( base, 2, e_adj.dart() );
                            if( maybe_phi2.has_value() )
                            {
                                const double barycentric = ( value - right_val ) / ( left_val - right_val );
                                mIntersectionPositions.emplace( Vertex( maybe_phi2.value() ), barycentric );
                            }
                        }
                        return base_dim > 2; // in 2d, there is only ever one edge to find
                    }
                }
                return true;
            } );
            return true;
        } );
    }
}

std::optional<Dart> LevelSetCMap::phi( const int i, const Dart& d ) const
{
    if( i > 1 ) return std::nullopt;
    const auto find_mark_in_phi1_chain = [this]( const Dart& a ) {
        for( Dart d_next = topology::phi( mUnderlyingMap, 1, a ).value(); d_next != a;
             d_next = topology::phi( mUnderlyingMap, 1, d_next ).value() )
        {
            if( mIntersectedEdges.isMarked( Edge( d_next ) ) ) return std::optional<Dart>( d_next );
        }
        throw std::runtime_error( "No next marked dart found!" );
    };

    if( i == 1 )
        return find_mark_in_phi1_chain( d ).and_then(
            [&]( const auto& phi1 ) { return topology::phi( mUnderlyingMap, 2, phi1 ); } );
    else
        return topology::phi( mUnderlyingMap, 2, d ).and_then( [&]( const auto& phi2 ) {
            return find_mark_in_phi1_chain( phi2 );
        } );
}

bool LevelSetCMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    if( mUnderlyingMap.dim() > 2 ) throw std::runtime_error( "Level sets on 3d maps not yet supported" );
    return iterateCellsWhile( 1, [&]( const Edge& e ) { return callback( e.dart() ); } );
}

bool LevelSetCMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( mUnderlyingMap.dim() > 2 ) throw std::runtime_error( "Level sets on 3d maps not yet supported" );
    if( cell_dim <= 1 )
    {
        for( const auto& pr : mIntersectionPositions )
        {
            if( not callback( Cell( pr.first.dart(), cell_dim ) ) ) return false;
        }
        return true;
    }
    else
    {
        return true;
    }
}

std::optional<IndexingFunc> LevelSetCMap::indexing( const uint cell_dim ) const
{
    return mUnderlyingMap.indexing( cell_dim + 1 )
        .and_then( [&]( const IndexingFunc underlying_func ) -> std::optional<IndexingFunc> {
            return [underlying_func]( const Cell& c ) { return underlying_func( Cell( c.dart(), c.dim() + 1 ) ); };
        } );
}