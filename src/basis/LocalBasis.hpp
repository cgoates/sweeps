#pragma once
#include <ParentDomain.hpp>

namespace basis
{
    enum class BasisType {
        Bernstein,
        BSpline
    };

    class BarycentricBasis
    {
        public:
        BasisType type = BasisType::Bernstein;
        size_t num_functions;
        //std::optional< std::vector<double> > mKnotVector;
    };

    size_t degree( const BarycentricBasis& bb );

    BarycentricBasis bernsteinBasis( const size_t degree );

    class LocalBasis
    {
        public:
        param::ParentDomain mParentDomain;
        SmallVector<BarycentricBasis, 3> mBasisGroups;
    };

    LocalBasis bernsteinSimplex( const size_t dim, const size_t degree );
    LocalBasis bernsteinCube( const size_t dim, const size_t degree );
}