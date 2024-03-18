#pragma once
#include <array>
#include <ostream>
#include <Eigen/Dense>

class VertexId
{
    public:
    using Type = uint64_t;
    VertexId( const Type& id );
    VertexId() {}
    Type id() const { return mId; }
    bool operator<( const VertexId& o ) const { return id() < o.id(); }
    bool operator>( const VertexId& o ) const { return id() > o.id(); }
    bool operator==( const VertexId& o ) const { return id() == o.id(); }
    friend std::ostream& operator<<( std::ostream& o, const VertexId& vid )
    {
        o << vid.mId;
        return o;
    }

    private:
    Type mId;
};

class Simplex
{
    public:
    Simplex( const VertexId& v0 );
    Simplex( const VertexId& v0, const VertexId& v1 );
    Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2 );
    Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2, const VertexId& v3 );
    template< typename It > requires std::input_iterator< It >
    Simplex( const It& v_list_begin, const size_t size ) :
        mDim( size - 1 ),
        mVertexIds( {0,0,0,0} )
    {
        if( size < 1 or size > 4 ) throw( "Bad Simplex size" );
        for( size_t i = 0; i < size; i++ ) mVertexIds.at( i ) = *std::next( v_list_begin, i );
    }

    const VertexId& vertex( const size_t n ) const;
    size_t dim() const { return mDim; }

    friend std::ostream& operator<<( std::ostream& os, const Simplex& s )
    {
        os << "{";
        for( size_t i = 0; i <= s.mDim; i++ )
        {
            os << s.mVertexIds.at( i ).id();
            if( i < s.mDim ) os << ", ";
        }
        os << "}";
        return os;
    }

    bool operator==( const Simplex& o ) const
    {
        if( dim() != o.dim() ) return false;
        for( size_t i = 0; i <= o.dim(); i++ )
        {
            if( vertex( i ) != o.vertex( i ) ) return false;
        }
        return true;
    }

    private:
    size_t mDim;
    std::array<VertexId, 4> mVertexIds;
};

class BarycentricPoint
{
    public:
    BarycentricPoint( const Simplex& s, const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3>& pt ) :
        simplex( s ),
        point( pt )
    {}

    Simplex simplex;
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3> point;

    // HACK
    bool operator==( const BarycentricPoint& o ) const
    {
        return o.simplex == simplex;
    }
};

BarycentricPoint addVertex( const BarycentricPoint lower_dim_pt, const VertexId& vid, const size_t new_point_pos );