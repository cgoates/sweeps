#include <BasisComplex1d.hpp>

namespace basis
{
    BasisComplex1d::BasisComplex1d( const param::ParametricAtlas1d& pa, const uint degree )
        : mAtlas( pa ), mParentBasis( bernsteinSimplex( 1, degree ) )
    {}

    const param::ParametricAtlas1d& BasisComplex1d::parametricAtlas() const
    {
        return mAtlas;
    }

    ParentBasis BasisComplex1d::parentBasis( const topology::Cell& ) const
    {
        return mParentBasis;
    }

    BasisComplex1d reduceDegree( const BasisComplex1d& bc )
    {
        return BasisComplex1d( bc.parametricAtlas(), bc.defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) );
    }
}