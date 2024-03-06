#include <CombinatorialMapMethods.hpp>
#include <queue>
#include <GlobalDartMarker.hpp>

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
                             const std::function<bool( const Dart& )>& callback )
    {
        // for each of the phi ops in the orbit, take that phi op and push it back into the queue.
        const int cell_dim = c.dim();
        constexpr int max_topo_dim = 3; // Maximum supported topology dimension
        std::queue<Dart> dart_queue;
        dart_queue.push( c.dart() );
        GlobalDartMarker m( map );

        const auto add_to_queue = [&]( const std::optional<Dart>& d ) {
            if( d.has_value() )
            {
                if( not m.isMarked( d.value() ) )
                {
                    dart_queue.push( d.value() );
                    m.mark( d.value() );
                }
            }
        };
        while( not dart_queue.empty() )
        {
            const Dart curr_d = dart_queue.front();
            dart_queue.pop();
            if( not callback( curr_d ) ) return false;
            if( cell_dim == 0 )
            {
                for( int j = 2; j <= max_topo_dim; j++ )
                {
                    add_to_queue( phi( map, { j, 1 }, curr_d ) );
                    add_to_queue( phi( map, { -1, j }, curr_d ) );
                }
            }
            else
            {
                for( int j = 1; j <= max_topo_dim; j++ )
                {
                    if( j != cell_dim ) add_to_queue( phi( map, j, curr_d ) );
                }
            }
        }
        return true;
    }

    size_t cellCount( const CombinatorialMap& map, const uint cell_dim )
    {
        size_t i = 0;
        map.iterateCellsWhile( cell_dim, [&]( const auto& ) {
            i++;
            return true;
        } );
        return i;
    }

} // namespace topology