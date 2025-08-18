#include <VectorConformingBasisComplex.hpp>

namespace basis
{
    VectorConformingBasisComplex::VectorConformingBasisComplex( const std::shared_ptr<const BasisComplex>& primal_complex,
                                                              const ConformingType type )
        : mPrimalComplex( primal_complex ), mConformingType( type )
    {}

    const param::ParametricAtlas& VectorConformingBasisComplex::parametricAtlas() const
    {
        return mPrimalComplex->parametricAtlas();
    }

    ParentBasis VectorConformingBasisComplex::parentBasis( const topology::Cell& c ) const
    {
        const ParentBasis primal = mPrimalComplex->parentBasis( c );
        const size_t param_dim = dim( primal.mParentDomain );
        if( not isCartesian( primal.mParentDomain ) or ( 2 != param_dim and 3 != param_dim ) )
            throw std::runtime_error( "VectorConformingBasisComplex only supports 2d or 3d cube-like cells" );
        
        SmallVector<size_t, 3> degrees;
        for( const auto& group : primal.mBasisGroups )
        {
            degrees.push_back( group.degrees.at( 0 ) );
        }

        if( mConformingType == ConformingType::Curl )
        {
            return ParentBasis{ primal.mParentDomain, { curlConformingBernsteinBasis( degrees ) } };
        }
        else if( mConformingType == ConformingType::Divergence )
        {
            return ParentBasis{ primal.mParentDomain, { divConformingBernsteinBasis( degrees ) } };
        }
        else
        {
            throw std::runtime_error( "Unknown conforming type" );
        }
    }
}