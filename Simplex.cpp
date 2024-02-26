#include <Simplex.hpp>
#include <cassert>

VertexId::VertexId( const VertexId::Type& id ) : mId( id ) {}

Simplex::Simplex( const VertexId& v0 )
    : mDim( 0 ), mVertexIds( { v0, 0, 0, 0 } )
{}

Simplex::Simplex( const VertexId& v0, const VertexId& v1 )
    : mDim( 1 ), mVertexIds( { v0, v1, 0, 0 } )
{}

Simplex::Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2 )
    : mDim( 2 ), mVertexIds( { v0, v1, v2, 0 } )
{}

Simplex::Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2, const VertexId& v3 )
    : mDim( 3 ), mVertexIds( { v0, v1, v2, v3 } )
{}

const VertexId& Simplex::vertex( const size_t n ) const
{
    assert( n <= mDim );
    return mVertexIds.at( n );
}

BarycentricPoint addVertex( const BarycentricPoint lower_dim_pt, const VertexId& vid, const size_t new_point_pos )
{
    if( new_point_pos > lower_dim_pt.simplex.dim() + 1 ) throw( "Bad new point position" );

    std::vector<VertexId> new_simplex;
    new_simplex.reserve( lower_dim_pt.simplex.dim() + 2 );

    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3> new_point( lower_dim_pt.simplex.dim() + 2 );

    for( size_t i = 0, j = 0; i <= lower_dim_pt.simplex.dim(); i++, j++ )
    {
        if( j == new_point_pos )
        {
            new_simplex.push_back( vid );
            new_point( j ) = 0;
            j++;
        }
        new_simplex.push_back( lower_dim_pt.simplex.vertex( i ) );
        new_point( j ) = lower_dim_pt.point( i );
    }
    return BarycentricPoint( Simplex( new_simplex.begin(), new_simplex.size() ), new_point );
}