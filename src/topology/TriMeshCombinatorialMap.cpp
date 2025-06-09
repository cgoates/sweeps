#include <Dart.hpp>
#include<TriMeshCombinatorialMap.hpp>
#include<queue>
#include<GlobalCellMarker.hpp>
#include <Logging.hpp>

namespace topology
{
constexpr int darts_per_tri = 3;

inline std::array<VertexId, 2> sortEdgeVerts( const VertexId& v1, const VertexId& v2 )
{
    if( v1 > v2 ) return std::array<VertexId, 2>{ v2, v1 };
    else return std::array<VertexId, 2>{ v1, v2 };
}

Dart TriMeshCombinatorialMap::dartOfTri( const uint tri_id ) const
{
    return Dart( tri_id * darts_per_tri ); 
}

Dart TriMeshCombinatorialMap::phi1_1( const int i, const Dart& d ) const
{
    const Dart::IndexType face_local_id = d.id() % darts_per_tri;
    return Dart( ( face_local_id + darts_per_tri + i ) % darts_per_tri
                 + d.id() - face_local_id );
}

TriMeshCombinatorialMap::TriMeshCombinatorialMap( const SimplicialComplex& complex ) :
    mSimplicialComplex( complex ),
    mPhi2s( maxDartId() + 1, maxDartId() + 1 ),
    mEdgeIds( maxDartId() + 1 ),
    mVertexDarts( complex.points.size() )
{
    const auto vertex_ids = indexingOrError( *this, 0 );

    std::map< std::array< VertexId, 2 >, std::pair< Dart, VertexId > > to_be_paired;

    size_t edge_ii = 0;
    for( size_t i = 0; i < complex.simplices.size(); i++ )
    {
        const Simplex& simplex = complex.simplices.at( i );
        if( simplex.dim() != 2 ) throw std::runtime_error( "Unsupported simplex dimension for TriMeshCombinatorialMap" );

        const Dart simplex_dart = dartOfTri( i );
        for( Dart::IndexType local_dart_id = 0; local_dart_id < darts_per_tri; local_dart_id++ )
        {
            const Dart d( simplex_dart.id() + local_dart_id );
            const VertexId dart_vertex = simplex.vertex( local_dart_id );
            const VertexId end_of_dart_vertex = simplex.vertex( ( local_dart_id + 1 ) % darts_per_tri );
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
            mVertexDarts.at( simplex.vertex( local_dart_id ).id() ) = d;
        }
    }
    mNumEdges = edge_ii;
}

Dart::IndexType TriMeshCombinatorialMap::maxDartId() const
{
    return mSimplicialComplex.simplices.size() * darts_per_tri - 1;
}

std::optional<Dart> TriMeshCombinatorialMap::phi( const int i, const Dart& d ) const
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


bool TriMeshCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    const Dart::IndexType max_dart_id = maxDartId();
    for( Dart::IndexType d = 0; d <= max_dart_id; d++ )
    {
        if( not callback( Dart( d ) ) ) return false;
    }
    return true;
}

bool TriMeshCombinatorialMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 2 )
    {
        for( size_t i = 0; i < mSimplicialComplex.simplices.size(); i++ )
        {
            if( not callback( Face( dartOfTri( i ) ) ) ) return false;
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

std::optional<IndexingFunc> TriMeshCombinatorialMap::indexing( const uint cell_dim ) const
{
    if( cell_dim == 0 )
    {
        return [this]( const Vertex& v ){
            const size_t simplex_id = v.dart().id() / darts_per_tri;
            const Simplex& simp = mSimplicialComplex.simplices.at( simplex_id );
            const VertexId::Type local_vert = v.dart().id() % darts_per_tri;
            return simp.vertex( local_vert ).id();
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
            return elem.dart().id() / darts_per_tri;
        };
    }
    return std::nullopt;
}

std::optional<size_t> TriMeshCombinatorialMap::cellCount( const uint cell_dim ) const
{
    switch( cell_dim )
    {
        case 0: return mSimplicialComplex.points.size();
        case 1: return mNumEdges;
        case 2: return mSimplicialComplex.simplices.size();
        default: return 0;
    }
}

}