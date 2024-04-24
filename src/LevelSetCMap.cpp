#include <LevelSetCMap.hpp>
#include <queue>
#include <Simplex.hpp>
#include <GlobalCellMarker.hpp>

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

    const auto underlying_vertex_ids = indexingOrError( base, 0 );

    const auto add_intersections = [&]( const Dart& d, const double coord ) {
        if( base_dim <= 2 ) mIntersectionPositions.emplace( Vertex( d ), coord );
        else
        {
            // This will actually work for 2d base as well, but the above is an optimization
            const size_t vertex_id = underlying_vertex_ids( Vertex( d ) );
            iterateDartsOfCell( base, Edge( d ), [&]( const Dart& other_d ) {
                const Vertex other_v( other_d );
                if( underlying_vertex_ids( other_v ) == vertex_id )
                    mIntersectionPositions.emplace( other_v, coord );
                return true;
            } );
        }
    };

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
                            add_intersections( e_adj.dart(), barycentric );
                        }
                        else
                        {
                            const auto maybe_phi2 = topology::phi( base, 2, e_adj.dart() );
                            if( maybe_phi2.has_value() )
                            {
                                const double barycentric = ( value - right_val ) / ( left_val - right_val );
                                add_intersections( maybe_phi2.value(), barycentric );
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
    else if( i == -1 )
        return topology::phi( mUnderlyingMap, 2, d ).and_then( [&]( const auto& phi2 ) {
            return find_mark_in_phi1_chain( phi2 );
        } );
    else if( i == 2 and dim() > 1 )
    {
        return topology::phi( mUnderlyingMap, 3, d ).and_then( [&]( const auto& phi3 ) {
            return find_mark_in_phi1_chain( phi3 );
        } );
    }
    else throw std::runtime_error( "Bad phi operation" );
}

bool LevelSetCMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    for( const auto& pr : mIntersectionPositions )
    {
        if( not callback( pr.first.dart() ) ) return false;
    }
    return true;
}

bool LevelSetCMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 0 and dim() == 1 )
    {
        for( const auto& pr : mIntersectionPositions )
        {
            if( not callback( pr.first ) ) return false;
        }
        return true;
    }
    else
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
}

std::optional<IndexingFunc> LevelSetCMap::indexing( const uint cell_dim ) const
{
    return mUnderlyingMap.indexing( cell_dim + 1 )
        .and_then( [&]( const IndexingFunc underlying_func ) -> std::optional<IndexingFunc> {
            return [underlying_func]( const Cell& c ) { return underlying_func( Cell( c.dart(), c.dim() + 1 ) ); };
        } );
}

double LevelSetCMap::intersectionPosition( const topology::Vertex& v ) const
{
    return mIntersectionPositions.at( v );
}

Cell LevelSetCMap::underlyingCell( const Cell& c ) const
{
    return Cell( c.dart(), c.dim() + 1 );
}

namespace topology
{
    std::function<Eigen::Vector3d( const Vertex& )> levelSetVertexPositions(
        const LevelSetCMap& level, const std::function<Eigen::Vector3d( const Vertex& )>& underlying_positions )
    {
        return [&, underlying_positions]( const Vertex& v ) -> Eigen::Vector3d {
            const double s = level.intersectionPosition( v );
            const Edge underlying_e = level.underlyingCell( v );
            const Vertex left_v( underlying_e.dart() );
            const Vertex right_v( phi( level.underlyingMap(), 1, underlying_e.dart() ).value() );
            return ( 1 - s ) * underlying_positions( left_v ) + s * underlying_positions( right_v );
        };
    }
}