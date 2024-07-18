#include <DivConfBasisComplex.hpp>

namespace basis
{
    DivConfBasisComplex::DivConfBasisComplex( const BasisComplex& primal_complex ) : mPrimalComplex( primal_complex ) {}

    const param::ParametricAtlas& DivConfBasisComplex::parametricAtlas() const
    {
        return mPrimalComplex.parametricAtlas();
    }

    ParentBasis DivConfBasisComplex::parentBasis( const topology::Cell& c ) const
    {
        const ParentBasis primal = mPrimalComplex.parentBasis( c );
        const size_t param_dim = dim( primal.mParentDomain );
        if( not isCartesian( primal.mParentDomain ) or ( 2 != param_dim and 3 != param_dim ) )
            throw std::runtime_error( "DivConfBasisComplex only supports 2d or 3d cube-like cells" );
        
        SmallVector<size_t, 3> degrees;
        for( const auto& group : primal.mBasisGroups )
        {
            degrees.push_back( group.degrees.at( 0 ) );
        }

        return ParentBasis{ primal.mParentDomain, { divConformingBernsteinBasis( degrees ) } };
    }
}