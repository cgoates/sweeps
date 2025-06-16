#include <MultiPatchDecomposition.hpp>
#include <MultiPatchCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <set>
#include <iostream>

namespace topology
{
    IndexingCellMarker multiPatchCorners( const CombinatorialMap& unstructured_cmap )
    {
        // Only used for markers; remove if necessary.
        const auto vert_indexing = indexingOrError( unstructured_cmap, 0 );
        IndexingCellMarker corners( vert_indexing, 0 );
        std::set<Vertex> extraordinary_vertices;
        iterateCellsWhile( unstructured_cmap, 0, [&]( const Vertex& v ) {
            if( isExtraordinary( unstructured_cmap, v ) )
            {
                extraordinary_vertices.insert( v );
                corners.mark( v );
            }
            return true;
        } );

        IndexingCellMarker trace_verts( vert_indexing, 0 );
        IndexingCellMarker trace_edges( indexingOrError( unstructured_cmap, 1 ), 1 );
        for( const auto& v : extraordinary_vertices )
        {
            iterateAdjacentCells( unstructured_cmap, v, 1, [&]( const Edge& e ) {
                if( onBoundary( unstructured_cmap, e.dart() ) ) return true;
                if( trace_edges.isMarked( e ) ) return true; // Already traced this edge.

                Dart curr_d = e.dart();
                if( vert_indexing( v ) != vert_indexing( Vertex( curr_d ) ) )
                {
                    const auto maybe_phi = phi( unstructured_cmap, 2, curr_d );
                    if( maybe_phi.has_value() )
                    {
                        curr_d = maybe_phi.value();
                    }
                    else
                    {
                        return true; // This is a boundary edge, so we don't need to trace it.
                    }
                }

                do
                {
                    trace_edges.mark( Edge( curr_d ) );
                    if( trace_verts.isMarked( Vertex( curr_d ) ) )
                        corners.mark( Vertex( curr_d ) );
                    else
                        trace_verts.mark( Vertex( curr_d ) );

                    const auto maybe_next_d = phi( unstructured_cmap, {1, 2, 1}, curr_d );
                    if( not maybe_next_d.has_value() )
                    {
                        // This is the boundary.
                        // Add the next vertex to corners.
                        corners.mark( Vertex( phi( unstructured_cmap, 1, curr_d ).value() ) );
                        break;
                    }
                    curr_d = maybe_next_d.value();
                } while( not isExtraordinary( unstructured_cmap, Vertex( curr_d ) ) );
                return true;
            } );
        }
        return corners;
    }

    std::pair<Dart, size_t> countInDirection( const CombinatorialMap& cmap, const Dart& d, const IndexingCellMarker& corners )
    {
        const std::vector<int> phi_ops = std::vector<int>{1, 2, 1};
        size_t count = 0;
        Dart curr_d = d;
        while( true )
        {
            ++count;
            const Dart next_v_dart = phi( cmap, 1, curr_d ).value();
            if( corners.isMarked( Vertex( next_v_dart ) ) )
            {
                return { next_v_dart, count };
            }
            const auto maybe_next_d = phi( cmap, phi_ops, curr_d );
            if( not maybe_next_d.has_value() )
            {
                throw std::runtime_error( "Reached the boundary but not a corner." );
            }

            curr_d = maybe_next_d.value();
        }
        return {d, count}; // This should never be reached.
    }

    MultiPatchDecomposition multiPatchDecomposition( const CombinatorialMap& unstructured_cmap )
    {
        const IndexingCellMarker corners = multiPatchCorners( unstructured_cmap );
        std::map<size_t, std::shared_ptr<const CombinatorialMap1d>> cmap_1ds;
        std::map<std::pair<size_t, size_t>, std::shared_ptr<const TPCombinatorialMap>> cmap_2ds;

        const auto find_or_add_1d = [&]( const size_t count ) -> std::shared_ptr<const CombinatorialMap1d> {
            auto it = cmap_1ds.find( count );
            if( it != cmap_1ds.end() )
            {
                return it->second;
            }
            else
            {
                auto new_cmap_1d = std::make_shared<const CombinatorialMap1d>( count );
                cmap_1ds.emplace( count, new_cmap_1d );
                return new_cmap_1d;
            }
        };
        const auto find_or_add_2d = [&]( const size_t count_s, const size_t count_t ) -> std::shared_ptr<const TPCombinatorialMap> {
            auto it = cmap_2ds.find( { count_s, count_t } );
            if( it != cmap_2ds.end() )
            {
                return it->second;
            }
            else
            {
                const auto cmap_1d_s = find_or_add_1d( count_s );
                const auto cmap_1d_t = find_or_add_1d( count_t );
                auto new_cmap_2d = std::make_shared<const TPCombinatorialMap>( cmap_1d_s, cmap_1d_t );
                cmap_2ds.emplace( std::pair<size_t, size_t>{count_s, count_t}, new_cmap_2d );
                return new_cmap_2d;
            }
        };

        const auto tp_corner_darts = []( const std::shared_ptr<const TPCombinatorialMap>& cmap, const size_t count_s, const size_t count_t ) -> std::array<Dart, 4> {
            const Dart corner0 = cmap->flatten( Dart( 0 ), Dart( 0 ), TPCombinatorialMap::TPDartPos::DartPos0 );
            const Dart corner1 = cmap->flatten( Dart( count_s - 1 ), Dart( 0 ), TPCombinatorialMap::TPDartPos::DartPos1 );
            const Dart corner2 = cmap->flatten( Dart( count_s - 1 ), Dart( count_t - 1 ), TPCombinatorialMap::TPDartPos::DartPos2 );
            const Dart corner3 = cmap->flatten( Dart( 0 ), Dart( count_t - 1 ), TPCombinatorialMap::TPDartPos::DartPos3 );

            return { corner0, corner1, corner2, corner3 };
        };


        GlobalDartMarker dart_marker( unstructured_cmap );

        std::vector<std::shared_ptr<const TPCombinatorialMap>> constituents;
        std::vector<std::array<Dart, 4>> constituent_corners;

        std::map<Dart, std::pair<size_t, Dart>> unstructured_to_constituent_darts;
        std::map<std::pair<size_t, Dart>, Dart> constituent_to_unstructured_darts;
        iterateCellsWhile( unstructured_cmap, 0, [&]( const Vertex& v ) {
            if( corners.isMarked( v ) )
            {
                iterateDartsOfCell( unstructured_cmap, v, [&]( const Dart& corner0 ) {
                    if( dart_marker.isMarked( corner0 ) )
                    {
                        // This corner has already been processed.
                        return true;
                    }
                    size_t count_s, count_t;
                    Dart corner1, corner2, corner3;

                    std::tie( corner1, count_s ) = countInDirection( unstructured_cmap, corner0, corners );
                    std::tie( corner2, count_t ) = countInDirection( unstructured_cmap, corner1, corners );
                    corner3 = countInDirection( unstructured_cmap, corner2, corners ).first;

                    dart_marker.mark( corner0 );
                    dart_marker.mark( corner1 );
                    dart_marker.mark( corner2 );
                    dart_marker.mark( corner3 );
                    constituent_corners.push_back( { corner0, corner1, corner2, corner3 } );

                    const size_t constituent_id = constituents.size();
                    constituents.push_back( find_or_add_2d( count_s, count_t ) );

                    const std::array<Dart, 4> tp_corners = tp_corner_darts( constituents.back(), count_s, count_t );
                    unstructured_to_constituent_darts.emplace( corner0, std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 0 ) } );
                    unstructured_to_constituent_darts.emplace( corner1, std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 1 ) } );
                    unstructured_to_constituent_darts.emplace( corner2, std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 2 ) } );
                    unstructured_to_constituent_darts.emplace( corner3, std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 3 ) } );
                    constituent_to_unstructured_darts.emplace( std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 0 ) }, corner0 );
                    constituent_to_unstructured_darts.emplace( std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 1 ) }, corner1 );
                    constituent_to_unstructured_darts.emplace( std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 2 ) }, corner2 );
                    constituent_to_unstructured_darts.emplace( std::pair<size_t, Dart>{ constituent_id, tp_corners.at( 3 ) }, corner3 );

                    return true;
                } );
            }
            return true;
        } );

        std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>> connections;

        std::set<std::pair<size_t, size_t>> connected_constituents;
        const auto add_connected_pair = [&]( const size_t id1, const size_t id2, const Dart& dart1, const Dart& dart2 ) {

            std::pair<size_t, size_t> pair = std::make_pair( std::min( id1, id2 ), std::max( id1, id2 ) );
            if( connected_constituents.find( pair ) == connected_constituents.end() )
            {
                connections.emplace( std::make_pair( id1, dart1 ), std::make_pair( id2, dart2 ) );
                connected_constituents.insert( pair );
            }
        };

        for( size_t constituent_ii = 0; constituent_ii < constituents.size(); ++constituent_ii )
        {
            const auto& constituent = constituents.at( constituent_ii );
            const auto tp_submaps = tensorProductComponentCMaps( *constituent );
            const std::array<Dart, 4> corners = tp_corner_darts( constituent, cellCount( *tp_submaps.at( 0 ), 1 ), cellCount( *tp_submaps.at( 1 ), 1 ) );
            for( size_t i = 0; i < 4; ++i )
            {
                const Dart& tp_corner_dart = corners.at( i );
                const Dart& unstructured_corner_dart = constituent_to_unstructured_darts.at( { constituent_ii, tp_corner_dart } );
                const auto maybe_phi2 = phi( unstructured_cmap, 2, unstructured_corner_dart );
                if( maybe_phi2.has_value() )
                {
                    const Dart& next_dart = maybe_phi2.value();
                    const auto [other_constituent_ii, other_tp_corner_dart] = unstructured_to_constituent_darts.at( phi( unstructured_cmap, 1, next_dart ).value() );
                    add_connected_pair( constituent_ii, other_constituent_ii, tp_corner_dart, phi( *constituents.at( other_constituent_ii ), -1, other_tp_corner_dart ).value() );
                }
            }
        }

        return { constituents, connections, constituent_to_unstructured_darts.at( { 0, Dart( 0 ) } ) };
    }
        
}