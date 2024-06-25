#include <TPBasisComplex.hpp>

namespace basis
{
    TPBasisComplex::TPBasisComplex( const param::TPParametricAtlas& pa,
                                    const BasisComplex& source_complex,
                                    const BasisComplex1d& line_complex )
        : mAtlas( pa ), mSourceComplex( source_complex ), mLineComplex( line_complex )
    {}

    const param::TPParametricAtlas& TPBasisComplex::parametricAtlas() const
    {
        return mAtlas;
    }

    const ParentBasis TPBasisComplex::parentBasis( const topology::Cell& c ) const
    {
        if( c.dim() != mAtlas.cmap().dim() ) throw std::runtime_error( "parentBasis only accepts cells of dimension cmap.dim()" );

        const auto [source_dart, line_dart, tp_op] = mAtlas.cmap().unflatten( c.dart() );
        const ParentBasis source_basis =
            mSourceComplex.parentBasis( topology::Cell( source_dart, mSourceComplex.parametricAtlas().cmap().dim() ) );
        const ParentBasis line_basis =
            mLineComplex.parentBasis( topology::Cell( line_dart, mLineComplex.parametricAtlas().cmap().dim() ) );

        return tensorProduct( source_basis, line_basis );
    }
}