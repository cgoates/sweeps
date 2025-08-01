#include <CombinatorialMapMethods.hpp>
#include <SmallQueue.hpp>
#include <set>
#include <GlobalDartMarker.hpp>
#include <GlobalCellMarker.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <Logging.hpp>
#include <queue>

namespace topology
{

    std::optional<Dart> phi( const CombinatorialMap& map, const int phi_op, const Dart& d )
    {
        return map.phi( phi_op, d );
    }

    std::optional<Dart> phi( const CombinatorialMap& map, const std::vector<int>& phi_ops, const Dart& d )
    {
        std::optional<Dart> out( d );
        for( const int& phi_op : phi_ops )
        {
            out = out.and_then( [&]( const Dart& a ) { return phi( map, phi_op, a ); } );
        }
        return out;
    }

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             DartMarker auto& m,
                             const std::function<bool( const Dart& )>& callback )
    {
        // for each of the phi ops in the orbit, take that phi op and push it back into the queue.
        const int cell_dim = c.dim();
        const int topo_dim = map.dim();

        GrowableQueue<Dart, 300> dart_queue;
        dart_queue.push( c.dart() );
        m.mark( c.dart() );

        const auto add_to_queue = [&m,&dart_queue]( const std::optional<Dart>& d ) {
            if( d.has_value() and not m.isMarked( d.value() ) )
            {
                dart_queue.push( d.value() );
                m.mark( d.value() );
            }
        };
        while( not dart_queue.empty() )
        {
            const Dart curr_d = dart_queue.pop();
            if( not callback( curr_d ) ) return false;
            if( cell_dim == 0 )
            {
                for( int j = 2; j <= topo_dim; j++ )
                {
                    add_to_queue( phi( map, { j, 1 }, curr_d ) );
                    add_to_queue( phi( map, { -1, j }, curr_d ) );
                }
            }
            else
            {
                for( int j = 1; j <= topo_dim; j++ )
                {
                    if( j != cell_dim ) add_to_queue( phi( map, j, curr_d ) );
                }
            }
        }
        return true;
    }
    template
    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             GlobalDartMarker& m,
                             const std::function<bool( const Dart& )>& callback );
    template
    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             LocalDartMarker& m,
                             const std::function<bool( const Dart& )>& callback );

    bool iterateDartsOfCell( const CombinatorialMap& map,
                             const Cell& c,
                             const std::function<bool( const Dart& )>& callback )
    {
        LocalDartMarker m;
        return iterateDartsOfCell( map, c, m, callback );
    }

    bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                       const Cell& c,
                                       const int restrict_dim,
                                       DartMarker auto& m,
                                       const std::function<bool( const Dart& )>& callback )
    {
        // for each of the phi ops in the orbit, take that phi op and push it back into the queue.
        const int cell_dim = c.dim();
        const int topo_dim = map.dim();

        SmallQueue<Dart, 300> dart_queue;
        dart_queue.push( c.dart() );
        m.mark( c.dart() );

        const auto add_to_queue = [&]( const std::optional<Dart>& d ) {
            if( d.has_value() and not m.isMarked( d.value() ) )
            {
                dart_queue.push( d.value() );
                m.mark( d.value() );
            }
        };
        while( not dart_queue.empty() )
        {
            const Dart curr_d = dart_queue.pop();
            if( not callback( curr_d ) ) return false;
            if( cell_dim == 0 )
            {
                if( restrict_dim == 1 ) continue;
                for( int j = 2; j <= topo_dim; j++ )
                {
                    if( j == restrict_dim ) continue;
                    add_to_queue( phi( map, { j, 1 }, curr_d ) );
                    add_to_queue( phi( map, { -1, j }, curr_d ) );
                }
            }
            else
            {
                for( int j = 1; j <= topo_dim; j++ )
                {
                    if( j != cell_dim and j != restrict_dim ) add_to_queue( phi( map, j, curr_d ) );
                }
            }
        }
        return true;
    }
    template bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                                const Cell& c,
                                                const int restrict_dim,
                                                GlobalDartMarker& m,
                                                const std::function<bool( const Dart& )>& callback );
    template bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                                const Cell& c,
                                                const int restrict_dim,
                                                LocalDartMarker& m,
                                                const std::function<bool( const Dart& )>& callback );

    bool iterateDartsOfRestrictedCell( const CombinatorialMap& map,
                                       const Cell& c,
                                       const int restrict_dim,
                                       const std::function<bool( const Dart& )>& callback )
    {
        LocalDartMarker m;
        return iterateDartsOfRestrictedCell( map, c, restrict_dim, m, callback );
    }

    bool iterateAdjacentCells( const CombinatorialMap& map,
                               const Cell& c,
                               const uint cell_dim,
                               const std::function<bool( const Cell& )>& callback )
    {
        // We need an extra operation for edges of vertices in map.dim() < 3
        const bool edges_of_vert_in_lowD = c.dim() == 0 and cell_dim == 1 and map.dim() <= 2;
        const bool verts_of_edge_in_lowD = c.dim() == 1 and cell_dim == 0 and map.dim() <= 2;

        const auto iterate = [&]( auto& m) {
            const auto mark_and_callback = [&]( const topology::Cell& cell ) {
                if( not m.isMarked( cell ) )
                {
                    m.mark( map, cell );
                    return callback( cell );
                }
                return true;
            };

            return iterateDartsOfCell( map, c, [&]( const Dart& d ){
                const Cell c_adj( d, cell_dim );
                if( not mark_and_callback( c_adj ) ) return false;
                if( edges_of_vert_in_lowD )
                {
                    const auto phi_1 = phi( map, -1, d );
                    if( phi_1.has_value() )
                    {
                        const topology::Edge c_adj_2( phi_1.value() );
                        if( not mark_and_callback( c_adj_2 ) ) return false;
                    }
                }
                else if( verts_of_edge_in_lowD )
                {
                    const auto phi1 = phi( map, 1, d );
                    if( phi1.has_value() )
                    {
                        const topology::Vertex c_adj_2( phi1.value() );
                        if( not mark_and_callback( c_adj_2 ) ) return false;
                    }
                }
                return true;
            } );
        };

        const auto indexing = map.indexing( cell_dim );
        if( indexing.has_value() )
        {
            IndexingCellMarker m( *indexing, cell_dim );
            return iterate( m );
        }
        else
        {
            LocalCellMarker m( cell_dim );
            return iterate( m );
        }
    }

    bool iterateAdjacentCellsOfRestrictedCell( const CombinatorialMap& map,
                                               const Cell& c,
                                               const int restrict_dim,
                                               const uint cell_dim,
                                               const std::function<bool( const Cell& )>& callback )
    {
        // We need an extra operation for edges of vertices in map.dim() < 3
        const bool edges_of_vert_in_lowD = c.dim() == 0 and cell_dim == 1 and map.dim() <= 2 and restrict_dim != 1;
        const bool verts_of_edge_in_lowD = c.dim() == 1 and cell_dim == 0 and map.dim() <= 2 and restrict_dim != 1;

        LocalCellMarker m( cell_dim );
        const auto mark_and_callback = [&]( const topology::Cell& cell ) {
            if( not m.isMarked( cell ) )
            {
                m.mark( map, cell );
                return callback( cell );
            }
            return true;
        };

        return iterateDartsOfRestrictedCell( map, c, restrict_dim, [&]( const Dart& d ){
            const Cell c_adj( d, cell_dim );
            if( not mark_and_callback( c_adj ) ) return false;
            if( edges_of_vert_in_lowD )
            {
                const auto phi_1 = phi( map, -1, d );
                if( phi_1.has_value() )
                {
                    const topology::Edge c_adj_2( phi_1.value() );
                    if( not mark_and_callback( c_adj_2 ) ) return false;
                }
            }
            else if( verts_of_edge_in_lowD )
            {
                const auto phi1 = phi( map, 1, d );
                if( phi1.has_value() )
                {
                    const topology::Vertex c_adj_2( phi1.value() );
                    if( not mark_and_callback( c_adj_2 ) ) return false;
                }
            }
            return true;
        } );
    }

    bool iterateCellsWhile( const CombinatorialMap& map,
                            const uint cell_dim,
                            const std::function<bool( const Cell& )>& callback )
    {
        return map.iterateCellsWhile( cell_dim, callback );
    }

    bool iterateDartsWhile( const CombinatorialMap& map,
                            const std::function<bool( const Dart& )>& callback )
    {
        return map.iterateDartsWhile( callback );
    }

    bool iterateReachableWhile( const CombinatorialMap& map,
                                const Dart& start_dart,
                                const std::vector<int>& phi_ops,
                                const std::function<bool( const Dart& )>& callback )
    {
        callback( start_dart );

        Dart curr_d = start_dart;
        while( true )
        {
            const auto next_d = phi( map, phi_ops, curr_d );
            if( not next_d.has_value() ) break;
            curr_d = next_d.value();
            if( curr_d == start_dart ) return true;
            if( not callback( curr_d ) ) return false;
        }

        std::vector<int> reverse_phi_ops( phi_ops.rbegin(), phi_ops.rend() );
        for( int& phi_op : reverse_phi_ops )
        {
            if( std::abs( phi_op ) == 1 )
            {
                phi_op = -phi_op;
            }
        }

        curr_d = start_dart;
        while( true )
        {
            const auto next_d = phi( map, reverse_phi_ops, curr_d );
            if( not next_d.has_value() ) break;
            curr_d = next_d.value();
            if( not callback( curr_d ) ) return false;
        }
        return true;
    }

    size_t cellCount( const CombinatorialMap& map, const uint cell_dim )
    {
        const auto cc = map.cellCount( cell_dim );
        if( cc.has_value() ) return cc.value();
        else
        {
            size_t i = 0;
            iterateCellsWhile( map, cell_dim, [&]( const auto& ) {
                i++;
                return true;
            } );
            return i;
        }
    }

    bool onBoundary( const CombinatorialMap& map, const Dart& d )
    {
        return not phi( map, map.dim(), d ).has_value();
    }

    bool boundaryAdjacent( const CombinatorialMap& map, const Cell& c )
    {
        return not iterateDartsOfCell( map, c, [&map]( const Dart& d ) {
            return not onBoundary( map, d );
        } );
    }

    std::optional<Dart> maybeBoundaryDart( const CombinatorialMap& map, const Cell& c )
    {
        std::optional<Dart> out;
        iterateDartsOfCell( map, c, [&]( const Dart& d ) {
            if( onBoundary( map, d ) )
            {
                out = d;
                return false;
            }
            return true;
        } );
        return out;
    }

    bool isExtraordinary( const CombinatorialMap& map, const Vertex& v )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "isExtraordinary only works for 2D maps" );
        size_t num_edges = 0;
        iterateAdjacentCells( map, v, 1, [&]( const auto& ) {
            num_edges++;
            return true;
        } );

        if( boundaryAdjacent( map, v ) )
        {
            return num_edges != 3;
        }
        else
        {
            return num_edges != 4;
        }
    }

    std::set<Dart> boundaryComponentDarts( const CombinatorialMap& map )
    {
        std::set<Dart> out;
        GlobalDartMarker m( map );
        iterateDartsWhile( map, [&]( const Dart& d ) {
            if( onBoundary( map, d ) and not m.isMarked( d ) )
            {
                out.insert( d );
                CombinatorialMapBoundary bdry( map, {d} );
                iterateDartsWhile( bdry, [&]( const Dart& bdry_d ) {
                    m.mark( bdry_d );
                    return true;
                } );
            }
            return true;
        } );
        return out;
    }

    Dart::IndexType lowestDartId( const CombinatorialMap& map, const Cell& c )
    {
        Dart out = c.dart();
        iterateDartsOfCell( map, c, [&]( const Dart& d ) {
            out = std::min( d, out );
            return true;
        } );
        return out.id();
    }

    void flood2d( const topology::CombinatorialMap& map,
                  const topology::Face& f,
                  const std::function<bool( const topology::Face& )>& stop_condition,
                  const std::function<void( const topology::Face& )>& mark_callback,
                  const std::function<void( const topology::Face& )>& callback )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Bad cmap dimension for flood2d" );
        std::queue<topology::Face> to_process;
        to_process.push( f );

        for( ; not to_process.empty(); to_process.pop() )
        {
            const topology::Face& curr_face = to_process.front();
            if( stop_condition( curr_face ) ) continue;
            callback( curr_face );
            mark_callback( curr_face );
            for( const auto& d : { phi( map, 2, curr_face.dart() ),
                                   phi( map, { 1, 2 }, curr_face.dart() ),
                                   phi( map, { -1, 2 }, curr_face.dart() ) } )
            {
                if( d.has_value() and not stop_condition( topology::Face( d.value() ) ) )
                {
                    to_process.push( topology::Face( d.value() ) );
                }
            }
        }
    }

    bool synchronizedFlood2d( const CombinatorialMap& map1,
                              const CombinatorialMap& map2,
                              const Dart& d1,
                              const Dart& d2,
                              const std::function<bool( const Dart&, const Dart& )>& callback )
    {
        if( map1.dim() != 2 or map2.dim() != 2 ) throw std::runtime_error( "Bad cmap dimension for synchronizedFlood2d" );
        std::queue<std::pair<Dart, Dart>> to_process;
        GlobalDartMarker m1( map1 );
        to_process.push( { d1, d2 } );
        m1.mark( d1 );

        for( ; not to_process.empty(); to_process.pop() )
        {
            const auto [curr_d1, curr_d2] = to_process.front();
            if( not callback( curr_d1, curr_d2 ) ) return false;
            for( size_t phi_op : { 1, 2 } )
            {
                const auto maybe_next_d1 = phi( map1, phi_op, curr_d1 );
                const auto maybe_next_d2 = phi( map2, phi_op, curr_d2 );
                if( maybe_next_d1.has_value() and maybe_next_d2.has_value() )
                {
                    const Dart next_d1 = maybe_next_d1.value();
                    const Dart next_d2 = maybe_next_d2.value();
                    if( not m1.isMarked( next_d1 ) )
                    {
                        m1.mark( next_d1 );
                        to_process.push( { next_d1, next_d2 } );
                    }
                }
                else if( maybe_next_d1.has_value() or maybe_next_d2.has_value() )
                {
                    throw std::runtime_error( "Synchronized flood 2d: provided cmaps do not have the same topology." );
                }
            }
        }
        return true;
    }

    IndexingFunc indexingOrError( const CombinatorialMap& map, const uint cell_dim )
    {
        return map.indexing( cell_dim ).or_else( [&]() -> std::optional<IndexingFunc> {
            throw std::runtime_error( "Indexing for dim " + std::to_string( cell_dim ) + " doesn't exist" );
        } ).value();
    }

    int eulerCharacteristic( const CombinatorialMap& map )
    {
        int out = 0;
        for( size_t i = 0; i <= map.dim(); i++ )
        {
            const int sign = i % 2 ? -1 : 1;
            out += sign * cellCount( map, i );
        }
        return out;
    }
} // namespace topology