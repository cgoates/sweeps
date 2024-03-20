#include <CombinatorialMapMethods.hpp>
#include <SmallQueue.hpp>
#include <set>
#include <GlobalDartMarker.hpp>
#include <GlobalCellMarker.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <Logging.hpp>

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
        // FIXME: This fails for edges of a vertex in 2d
        LocalCellMarker m( cell_dim );
        return iterateDartsOfCell( map, c, [&]( const Dart& d ){
            const Cell c_adj( d, cell_dim );
            if( not m.isMarked( c_adj ) )
            {
                m.mark( map, c_adj );
                if( not callback( c_adj ) ) return false;
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
} // namespace topology