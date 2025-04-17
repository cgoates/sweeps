#include <HierarchicalTPBasisComplex.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CustomEigen.hpp>

using namespace basis;

HierarchicalTPBasisComplex::HierarchicalTPBasisComplex(
    const std::shared_ptr<const param::HierarchicalTPParametricAtlas>& pa,
    const std::vector<std::shared_ptr<const TPBasisComplex>>& refinement_levels )
    : mAtlas( pa ), mRefinementLevels( refinement_levels )
{}

const param::HierarchicalTPParametricAtlas& HierarchicalTPBasisComplex::parametricAtlas() const
{
    return *mAtlas;
}

ParentBasis HierarchicalTPBasisComplex::parentBasis( const topology::Cell& c ) const
{
    if( c.dim() != mAtlas->cmap().dim() )
        throw std::invalid_argument( "HierarchicalTPBasisComplex::parentBasis only works for elements." );

    const auto [ level, level_d ] = mAtlas->cmap().unrefinedAncestorDartOfCell( c );
    return mRefinementLevels.at( level )->parentBasis( topology::Cell( level_d, c.dim() ) );
}