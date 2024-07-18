#pragma once
#include <ParentDomain.hpp>

namespace basis
{
    enum class BasisType {
        Bernstein,
        DivConformingBernstein,
        //BSpline
    };

    class BarycentricBasis
    {
        public:
        BasisType type = BasisType::Bernstein;
        SmallVector<size_t, 3> degrees; // Will be of size 1 except for div conforming. Size indicates number of parametric groups in basis group.
        //std::optional< std::vector<double> > mKnotVector;
    };

    size_t numFunctions( const BarycentricBasis& bb );

    BarycentricBasis bernsteinBasis( const size_t degree );
    BarycentricBasis divConformingBernsteinBasis( const size_t dim, const size_t degree );
    BarycentricBasis divConformingBernsteinBasis( const SmallVector<size_t, 3>& degrees );

    class ParentBasis
    {
        public:
        param::ParentDomain mParentDomain;
        SmallVector<BarycentricBasis, 3> mBasisGroups;
    };

    size_t numFunctions( const ParentBasis& pb );
    SmallVector<size_t, 3> degrees( const ParentBasis& pb );
    size_t numVectorComponents( const ParentBasis& pb );

    ParentBasis bernsteinSimplex( const size_t dim, const size_t degree );
    ParentBasis bernsteinCube( const size_t dim, const size_t degree );
    ParentBasis divConformingBernsteinCube( const size_t dim, const size_t primal_degree );
    ParentBasis tensorProduct( const ParentBasis& pb1, const ParentBasis& pb2 );
}