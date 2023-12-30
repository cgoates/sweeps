#pragma once
#include<array>
#include<vector>
#include<Eigen/Dense>

class VertexId
{
    public:
    using Type = uint64_t;
    VertexId( const Type& id );
    Type id() const { return mId; }
    bool operator<( const VertexId& o ) const
    {
        return id() < o.id();
    }
    friend std::ostream& operator<<( std::ostream& o, const VertexId& vid )
    {
        o << vid.mId;
        return o;
    }
    private:
    const Type mId;
};

// TODO: I don't love the constructors here.  It would be nice to be able to use initializer lists.
class Simplex
{
    public:
    Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2 );
    Simplex( const VertexId& v0, const VertexId& v1, const VertexId& v2, const VertexId& v3 );

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

    private:
    size_t mDim;
    std::array<VertexId, 4> mVertexIds;
};

class SimplicialComplex
{
    public:
    SimplicialComplex( const std::vector<Simplex>& simplices ) : mSimplices( simplices ) {}

    const Simplex& simplex( const size_t simplex_id ) const;
    size_t numSimplices() const { return mSimplices.size(); }
    const std::vector<Simplex>& simplices() const { return mSimplices; }
    private:
    std::vector<Simplex> mSimplices;
};

struct LinearField
{
    LinearField( const SimplicialComplex& complex, const Eigen::MatrixXd& field ) :
        mComplex( complex ),
        mVertexFieldValues( field )
    {}

    const SimplicialComplex& mComplex;
    Eigen::MatrixXd mVertexFieldValues;
};