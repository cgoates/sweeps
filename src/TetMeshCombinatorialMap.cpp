#include <Dart.hpp>
#include<TetMeshCombinatorialMap.hpp>
#include<queue>
#include<GlobalCellMarker.hpp>
#include <Logging.hpp>

namespace topology
{

// TODO: Put this in a namespace so it is file local only
constexpr int darts_per_tri = 3;
constexpr int darts_per_tet = 12;

std::array<VertexId, 3> sortTriVerts( const VertexId& v1, const VertexId& v2, const VertexId& v3 )
{
    std::array<VertexId, 3> out( { v1, v2, v3 } );
    if( out[0] > out[1] ) std::swap( out[0], out[1] );
    if( out[0] > out[2] ) std::swap( out[0], out[2] );
    if( out[1] > out[2] ) std::swap( out[1], out[2] );
    return out;
}

Dart TetMeshCombinatorialMap::dartOfTet( const uint tet_id ) const
{
    return Dart( tet_id * darts_per_tet ); 
}

inline Dart phi1_1( const int i, const Dart& d )
{
    const Dart::IndexType face_local_id = d.id() % darts_per_tri;
    return Dart( ( face_local_id + darts_per_tri + i ) % darts_per_tri
                 + d.id() - face_local_id );
}

TetMeshCombinatorialMap::TetMeshCombinatorialMap( const SimplicialComplex& complex ) :
    mSimplicialComplex( complex ),
    mPhi3s( ( maxDartId() + 1 ) / darts_per_tri, maxDartId() + 1 ),
    mFaceIds( maxDartId() + 1 ),
    mEdgeIds( maxDartId() + 1 ),
    mVertexDarts( complex.points.size() )
{
    // Loop through all tets
    // Loop through four faces of tet

    const auto make_connection = [&]( const VertexId& v_at_end_of_d, const Dart& d_on_f_opp ) {
        Dart phi3 = d_on_f_opp;
        while( vertexId( phi3 ) != v_at_end_of_d )
        {
            phi3 = phi1_1( 1, phi3 );
        }
        return phi3;
    };

    std::map< std::array< VertexId, 3 >, std::pair< Dart, VertexId > > to_be_paired;

    size_t face_ii = 0;
    for( size_t i = 0; i < complex.simplices.size(); i++ )
    {
        const Simplex& simplex = complex.simplices.at( i );
        if( simplex.dim() != 3 ) throw std::runtime_error( "Unsupported simplex dimension for TetMeshCombinatorialMap" );
        for( const auto& local_dart_id : mLocalFaceDarts )
        {
            const Dart d( dartOfTet( i ).id() + local_dart_id );
            const VertexId dart_vertex = simplex.vertex( mLocalVertices.at( local_dart_id ) );
            const VertexId end_of_dart_vertex = simplex.vertex( mLocalVertices.at( local_dart_id + 1 ) );
            const std::array<VertexId, 3> sorted_tri = sortTriVerts( dart_vertex,
                                                                     end_of_dart_vertex,
                                                                     simplex.vertex( mLocalVertices.at( local_dart_id + 2 ) ) );
            
            const auto [it, inserted] = to_be_paired.try_emplace( sorted_tri, d, end_of_dart_vertex );
            size_t face_id;
            if( inserted )
            {
                face_id = face_ii++;
            }
            else
            {
                face_id = mFaceIds.at( it->second.first.id() );
                // do the magic to connect two phis
                mPhi3s.at( d.id() / darts_per_tri ) = make_connection( end_of_dart_vertex, it->second.first );
                mPhi3s.at( it->second.first.id() / darts_per_tri ) = make_connection( it->second.second, d );
            }

            mFaceIds.at( d.id() ) = face_id;
            mFaceIds.at( d.id() + 1 ) = face_id;
            mFaceIds.at( d.id() + 2 ) = face_id;
        }
        for( size_t j = 0; j <= simplex.dim(); j++ )
        {
            mVertexDarts.at( simplex.vertex( j ).id() ) = dartOfTet( i ).id() + mLocalVertexDarts.at( j );
        }
    }
    mNumFaces = face_ii;


    // Build edge indexing.
    size_t edge_ii = 0;
    // Using iterateDartsWhile here instead of iterateCellsWhile for better performance.
    GlobalDartMarker m( *this );
    iterateDartsWhile( [&]( const Dart& d ){
        if( not m.isMarked( d ) )
        {
            const Edge e( d );
            iterateDartsOfCell( *this, e, m, [&]( const Dart& d ) {
                mEdgeIds.at( d.id() ) = edge_ii;
                return true;
            } );
            edge_ii++;
        }
        return true;
    } );
    mNumEdges = edge_ii;
}

Dart::IndexType TetMeshCombinatorialMap::maxDartId() const
{
    return mSimplicialComplex.simplices.size() * darts_per_tet - 1;
}

std::optional<Dart> TetMeshCombinatorialMap::phi( const int i, const Dart& d ) const
{
    switch( i )
    {
        case -1:
        case 1: return phi1_1( i, d ); break;
        case 2:
        {
            const Dart::IndexType tet_local_id = d.id() % darts_per_tet;
            return Dart( mPhi2s.at( tet_local_id ) + d.id() - tet_local_id );
        }
        break;
        case 3:
        {
            // The key dart on each face is the lowest dart number on the face.
            const Dart::IndexType face_local_id = d.id() % darts_per_tri;

            Dart out = mPhi3s.at( d.id() / darts_per_tri );
            if( out.id() > maxDartId() ) return {};

            for( Dart::IndexType new_local_id = 0; new_local_id < face_local_id; new_local_id++ )
            {
                out = phi1_1( -1, out );
            }
            return out;
        }
        default:
            throw std::invalid_argument( "Invalid phi operation for tet mesh combinatorial map" );
    }
}


bool TetMeshCombinatorialMap::iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const
{
    const Dart::IndexType max_dart_id = maxDartId();
    for( Dart::IndexType d = 0; d <= max_dart_id; d++ )
    {
        if( not callback( Dart( d ) ) ) return false;
    }
    return true;
}

bool TetMeshCombinatorialMap::iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const
{
    if( cell_dim == 3 )
    {
        for( size_t i = 0; i < mSimplicialComplex.simplices.size(); i++ )
        {
            if( not callback( Volume( dartOfTet( i ) ) ) ) return false;
        }
    }
    else if( cell_dim == 0 )
    {
        for( const Dart& d : mVertexDarts )
        {
            if( not callback( Vertex( d ) ) ) return false;
        }
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
    return true;
}

size_t TetMeshCombinatorialMap::elementId( const Volume& c ) const
{
    return c.dart().id() / darts_per_tet;
}

VertexId TetMeshCombinatorialMap::vertexId( const Vertex& c ) const
{
    const size_t simplex_id = elementId( Cell( c.dart(), 3 ) );
    const Simplex& simp = mSimplicialComplex.simplices.at( simplex_id );
    const VertexId::Type local_vert = mLocalVertices.at( c.dart().id() % darts_per_tet );
    return simp.vertex( local_vert );
}

size_t TetMeshCombinatorialMap::faceId( const Face& f ) const
{
    return mFaceIds.at( f.dart().id() );
}

size_t TetMeshCombinatorialMap::edgeId( const Edge& e ) const
{
    return mEdgeIds.at( e.dart().id() );
}

std::optional<size_t> TetMeshCombinatorialMap::cellCount( const uint cell_dim ) const
{
    switch( cell_dim )
    {
        case 0: return mSimplicialComplex.points.size();
        case 1: return mNumEdges;
        case 2: return mNumFaces;
        case 3: return mSimplicialComplex.simplices.size();
        default: throw std::runtime_error( "bad cell dim" );
    }
}

}