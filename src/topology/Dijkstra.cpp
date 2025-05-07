#include <Dijkstra.hpp>
#include <CombinatorialMapMethods.hpp>
#include <set>
#include <SimplexUtilities.hpp>

namespace topology
{
    // TODO: This should be replaced with a fibonacci heap for performance, or a binary heap.
    // See https://stackoverflow.com/questions/649640/
    template<class T, class Comp>
    class PriorityQueue
    {
        public:
        PriorityQueue( const Comp& c ) : mSet( c ) {}
        void push( const T& t, const double d )
        {
            mSet.insert( {d, t} );
        }

        T pop()
        {
            if( mSet.size() == 0 ) throw std::runtime_error( "Empty queue" );
            const auto front = *mSet.begin();
            mSet.erase( mSet.begin() );
            return front.second;
        }

        bool empty() const
        {
            return mSet.empty();
        }

        void updatePriority( const T& t, const double d )
        {
            std::erase_if( mSet, [&]( const auto& a ) { return a.second == t; } );
            mSet.insert( {d, t} );
        }

        std::set<std::pair<double, T>, Comp> mSet;
    };

    std::vector<Edge> shortestPath( const CombinatorialMap& map,
                                    const VertexPositionsFunc& vert_positions,
                                    const Vertex& start_vertex,
                                    const std::function<bool( const Vertex& )>& stop_condition,
                                    const bool interior_only )
    {
        return shortestPath(
            map,
            [&]( const Edge& e ) { return edgeLength( map, vert_positions, e ); },
            start_vertex,
            stop_condition,
            interior_only );
    }

    std::vector<Edge> shortestPath( const CombinatorialMap& map,
                                    const std::function<double( const Edge& )>& edge_lengths,
                                    const Vertex& start_vertex,
                                    const std::function<bool( const Vertex& )>& stop_condition,
                                    const bool interior_only )
    {
        const auto vert_ids = indexingOrError( map, 0 );
        // NOTE: Lots of extra storage here for topologies with non-contiguous indexing. maybe FIXME?
        const size_t num_verts = [&]() {
            size_t out = 0;
            iterateCellsWhile( map, 0, [&]( const auto& v ) {
                out = std::max( out, vert_ids( v ) );
                return true;
            } );
            return out + 1;
        }();
        std::vector<double> distances( num_verts, std::numeric_limits<double>::max() );
        std::vector<size_t> path_lengths( num_verts, 0 );
        std::vector<std::optional<Vertex>> prev( num_verts, std::nullopt );

        const auto comp = [&]( const std::pair<double, Vertex>& v1, const std::pair<double, Vertex>& v2 ) {
            if( v1.first == v2.first ) return vert_ids( v1.second ) < vert_ids( v2.second );
            else return v1.first < v2.first;
        };

        PriorityQueue<Vertex, decltype( comp )> Q( comp );

        distances.at( vert_ids( start_vertex ) ) = 0;
        path_lengths.at( vert_ids( start_vertex ) ) = 0;
        Q.push( start_vertex, 0.0 );

        std::optional<Vertex> finish_vert;

        while( not Q.empty() and not finish_vert )
        {
            const Vertex& v = Q.pop();
            if( stop_condition( v ) )
            {
                finish_vert.emplace( v );
                break;
            }
            iterateDartsOfCell( map, v, [&]( const Dart& d ) {
                if( interior_only and onBoundary( map, d ) ) return true;
                const Vertex v_next( phi( map, 1, d ).value() );
                const double new_dist = distances.at( vert_ids( v ) ) + edge_lengths( d );
                if( new_dist < distances.at( vert_ids( v_next ) ) )
                {
                    path_lengths.at( vert_ids( v_next ) ) = path_lengths.at( vert_ids( v ) ) + 1;
                    prev.at( vert_ids( v_next ) ).emplace( d );
                    distances.at( vert_ids( v_next ) ) = new_dist;
                    Q.updatePriority( v_next, new_dist );
                }
                return true;
            } );
        }

        if( not finish_vert )
            throw std::runtime_error( "Unable to find a path" );

        // Process the shortest path
        std::vector<Edge> path( path_lengths.at( vert_ids( finish_vert.value() ) ) );
        Vertex curr_v = finish_vert.value();
        for( size_t i = 0; i < path.size(); i++ )
        {
            curr_v = prev.at( vert_ids( curr_v ) ).value();
            path.at( path.size() - i - 1 ) = Edge( curr_v.dart() );
        }

        return path;
    }
}