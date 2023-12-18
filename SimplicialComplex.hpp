#pragma once
#include<array>
#include<vector>
#include<Eigen/Dense>

class VertexId
{
    public:
    using Type = uint64_t;
    VertexId( const Type& id );
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

    private:
    size_t mDim;
    std::array<VertexId, 4> mVertexIds;
};

class SimplicialComplex
{
    public:
    SimplicialComplex( const std::vector<Simplex>& simplices ) : mSimplices( simplices ) {}

    const Simplex& simplex( const size_t simplex_id ) const;
    private:
    std::vector<Simplex> mSimplices;
};

class LinearField
{
    public:
    LinearField( const SimplicialComplex& complex, const Eigen::MatrixXd& field ) :
        mComplex( complex ),
        mVertexFieldValues( field )
    {}

    private:
    const SimplicialComplex& mComplex;

    Eigen::MatrixXd mVertexFieldValues;
};