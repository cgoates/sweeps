#include <Dart.hpp>
#include <QuadMeshCombinatorialMap.hpp>
#include <queue>
#include <GlobalCellMarker.hpp>
#include <Logging.hpp>

namespace topology
{
constexpr int darts_per_quad = 4;

std::array<VertexId, 2> sortEdgeVerts( const VertexId& v1, const VertexId& v2 )
{
    if( v1 > v2 ) return std::array<VertexId, 2>{ v2, v1 };
    else return std::array<VertexId, 2>{ v1, v2 };
}

Dart QuadMeshCombinatorialMap::dartOfQuad( const uint quad_id ) const
{
    return Dart( quad_id * darts_per_quad ); 
}

Dart QuadMeshCombinatorialMap::phi1_1( const int i, const Dart& d ) const
{
    const Dart::IndexType face_local_id = d.id() % darts_per_quad;
    return Dart( ( face_local_id + darts_per_quad + i ) % darts_per_quad
                 + d.id() - face_local_id );
}

QuadMeshCombinatorialMap::QuadMeshCombinatorialMap( const std::vector<std::array<VertexId, 4>>& quads, const size_t num_vertices ) :
    mVerticesOfQuads( quads ),
    mPhi2s( maxDartId() + 1, maxDartId() + 1 ),
    mEdgeIds( maxDartId() + 1 ),
    mVertexDarts( num_vertices )
{
    const auto vertex_ids = indexingOrError( *this, 0 );

    std::map< std::array< VertexId, 2 >, std::pair< Dart, VertexId > > to_be_paired;

    size_t edge_ii = 0;
    for( size_t i = 0; i < quads.size(); i++ )
    {
        const std::array<VertexId, 4>& quad = quads.at( i );

        const Dart quad_dart = dartOfQuad( i );
        for( Dart::IndexType local_dart_id = 0; local_dart_id < darts_per_quad; local_dart_id++ )
        {
            const Dart d( quad_dart.id() + local_dart_id );
            const VertexId dart_vertex = quad.at( local_dart_id );
            const VertexId end_of_dart_vertex = quad.at( ( local_dart_id + 1 ) % darts_per_quad );
            const std::array<VertexId, 2> sorted_edge = sortEdgeVerts( dart_vertex, end_of_dart_vertex );

            const auto [it, inserted] = to_be_paired.try_emplace( sorted_edge, d, end_of_dart_vertex );
            size_t edge_id;
            if( inserted )
            {
                edge_id = edge_ii++;
            }
            else
            {
                edge_id = mEdgeIds.at( it->second.first.id() );
                // do the magic to connect two phis
                mPhi2s.at( d.id() ) = it->second.first;
                mPhi2s.at( it->second.first.id() ) = d;
            }

            mEdgeIds.at( d.id() ) = edge_id;
            mVertexDarts.at( quad.at( local_dart_id ).id() ) = d;
        }
    }
    mNumEdges = edge_ii;
}

Dart::IndexType QuadMeshCombinatorialMap::maxDartId() const
{
    return mVerticesOfQuads.size() * darts_per_quad - 1;
}

std::optional<Dart> QuadMeshCombinatorialMap::phi( const int i, const Dart& d ) const
{
    switch( i )
    {
        case -1:
        case 1: return phi1_1( i, d ); break;
        case 2:
        {
            const auto maybe_phi2 = mPhi2s.at( d.id() );
            if( maybe_phi2.id() > maxDartId() ) return std::nullopt;
            return maybe_phi2;
            break;
        }
        default:
            throw std::invalid_argument( "Invalid phi operation for tet mesh combinatorial map" );
    }
}


bool QuadMeshCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    const Dart::IndexType max_dart_id = maxDartId();
    for( Dart::IndexType d = 0; d <= max_dart_id; d++ )
    {
        if( not callback( Dart( d ) ) ) return false;
    }
    return true;
}

bool QuadMeshCombinatorialMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 2 )
    {
        for( size_t i = 0; i < mVerticesOfQuads.size(); i++ )
        {
            if( not callback( Face( dartOfQuad( i ) ) ) ) return false;
        }
    }
    else if( cell_dim == 0 )
    {
        for( const Dart& d : mVertexDarts )
        {
            if( not callback( Vertex( d ) ) ) return false;
        }
    }
    else if( cell_dim == 1 )
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
    return true;
}

std::optional<IndexingFunc> QuadMeshCombinatorialMap::indexing( const uint cell_dim ) const
{
    if( cell_dim == 0 )
    {
        return [this]( const Vertex& v ){
            const size_t quad_id = v.dart().id() / darts_per_quad;
            const std::array<VertexId, 4>& quad = mVerticesOfQuads.at( quad_id );
            const VertexId::Type local_vert = v.dart().id() % darts_per_quad;
            return quad.at( local_vert ).id();
        };
    }
    else if( cell_dim == 1 )
    {
        return [this]( const Edge& e ){
            return mEdgeIds.at( e.dart().id() );
        };
    }
    else if( cell_dim == 2 )
    {
        return []( const Face& elem ) {
            return elem.dart().id() / darts_per_quad;
        };
    }
    return std::nullopt;
}

std::optional<size_t> QuadMeshCombinatorialMap::cellCount( const uint cell_dim ) const
{
    switch( cell_dim )
    {
        case 0: return mVertexDarts.size();
        case 1: return mNumEdges;
        case 2: return mVerticesOfQuads.size();
        default: return 0;
    }
}

}