#include <SimplicialComplex.hpp>


VertexId::VertexId( const VertexId::Type& id ) : mId( id ) {}

Simplex::Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2 ) :
    mDim( 2 ),
    mVertexIds( { v0, v1, v2, 0 } )
{}

Simplex::Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2, const VertexId& v3 ) :
    mDim( 3 ),
    mVertexIds( { v0, v1, v2, v3 } )
{}

const VertexId& Simplex::vertex( const size_t n ) const
{
    assert( n <= mDim );
    return mVertexIds.at( n );
}

const Simplex& SimplicialComplex::simplex( const size_t simplex_id ) const
{
    return mSimplices.at( simplex_id );
}