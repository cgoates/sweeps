#include<TetMeshCombinatorialMap.hpp>
#include<cassert>

// TODO: Put this in a namespace so it is file local only
constexpr int darts_per_tri = 3;
constexpr int darts_per_tet = 12;

TetMeshCombinatorialMap::TetMeshCombinatorialMap( const SimplicialComplex& complex ) :
    mSimplicialComplex( complex )
{
    struct FaceCellPair
    {

    };

    // Loop through all tets
    // Loop through four faces of tet

    std::map< std::array< VertexId, 3 >, std::pair< std::array<VertexId, 3>, Dart > > to_be_paired;

    

}

inline Dart phi1_1( const int i, const Dart& d )
{
    const Dart::IndexType face_local_id = d.index() % darts_per_tri;
    return Dart( ( face_local_id + darts_per_tri + i ) % darts_per_tri
                 + d.index() - face_local_id );
}

std::optional<Dart> TetMeshCombinatorialMap::phi( const int i, const Dart& d ) const
{
    switch( i )
    {
        case -1:
        case 1: return phi1_1( i, d ); break;
        case 2:
        {
            const Dart::IndexType tet_local_id = d.index() % darts_per_tet;
            return Dart( mPhi2s.at( tet_local_id ) + d.index() - tet_local_id );
        }
        break;
        case 3:
        {
            // The key dart on each face is the lowest dart number on the face.
            const Dart::IndexType face_local_id = d.index() % darts_per_tri;
            const Dart key_dart = d.index() - face_local_id;

            const auto it = mPhi3s.find( key_dart );
            if( it == mPhi3s.end() ) return {};

            Dart out = it->second;
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

size_t TetMeshCombinatorialMap::elementId( const Cell& c ) const
{
    assert( c.dim() == 3 );
    return c.dart().index() / darts_per_tet;
}

VertexId TetMeshCombinatorialMap::vertexId( const Cell& c ) const
{
    assert( c.dim() == 0 );
    const size_t simplex_id = elementId( Cell( c.dart(), 3 ) );
    const Simplex& simp = mSimplicialComplex.simplex( simplex_id );
    const VertexId::Type local_vert = mLocalVertices.at( c.dart().index() % darts_per_tet );
    return simp.vertex( local_vert );
}