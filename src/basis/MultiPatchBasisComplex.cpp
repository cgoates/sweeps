#include <MultiPatchBasisComplex.hpp>

namespace basis
{
    MultiPatchBasisComplex::MultiPatchBasisComplex(
        const std::shared_ptr<const param::MultiPatchParametricAtlas>& p_atlas,
        const std::vector<std::shared_ptr<const TPBasisComplex>>& constituents )
        : mAtlas( p_atlas ), mSubComplexes( constituents )
    {}

    const param::MultiPatchParametricAtlas& MultiPatchBasisComplex::parametricAtlas() const
    {
        return *mAtlas;
    }

    ParentBasis MultiPatchBasisComplex::parentBasis( const topology::Cell& c ) const
    {
        const auto [patch_id, local_d] = mAtlas->cmap().toLocalDart( c.dart() );
        const topology::Cell local_c( local_d, c.dim() );
        return mSubComplexes.at( patch_id )->parentBasis( c );
    }
}