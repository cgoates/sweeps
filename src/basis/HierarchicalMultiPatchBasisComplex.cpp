#include <HierarchicalMultiPatchBasisComplex.hpp>

using namespace basis;

std::vector<std::shared_ptr<const HierarchicalTPBasisComplex>>
    initializeConstituents( const param::HierarchicalMultiPatchParametricAtlas& pa,
                            const std::vector<std::shared_ptr<const MultiPatchBasisComplex>>& refinement_levels )
{
    if( pa.refinementLevels().size() != refinement_levels.size() )
        throw std::invalid_argument( "Inconsistent refinement levels and cmap provided to HierarchicalMultiPatchBasisComplex." );

    std::vector<std::vector<std::shared_ptr<const TPBasisComplex>>> constituent_levels( pa.constituents().size() );
    for( const auto& level : refinement_levels )
    {
        if( level->constituents().size() != pa.constituents().size() )
            throw std::invalid_argument( "All refinement levels must have the same number of constituents." );
        for( size_t patch_ii = 0; patch_ii < pa.constituents().size(); patch_ii++ )
        {
            constituent_levels.at( patch_ii ).push_back( level->constituents().at( patch_ii ) );
        }
    }

    std::vector<std::shared_ptr<const HierarchicalTPBasisComplex>> out;
    out.reserve( pa.constituents().size() );
    for( size_t i = 0; i < pa.constituents().size(); i++ )
    {
        out.push_back( std::make_shared<HierarchicalTPBasisComplex>( pa.constituents().at( i ), constituent_levels.at( i ) ) );
    }
    return out;
}

HierarchicalMultiPatchBasisComplex::HierarchicalMultiPatchBasisComplex( const std::shared_ptr<const param::HierarchicalMultiPatchParametricAtlas>& pa,
                                    const std::vector<std::shared_ptr<const MultiPatchBasisComplex>>& refinement_levels ) :
    mAtlas( pa ), mRefinementLevels( refinement_levels ), mConstituents( initializeConstituents( *pa, refinement_levels ) )
{}

ParentBasis HierarchicalMultiPatchBasisComplex::parentBasis( const topology::Cell& c ) const
{
    const auto [patch_ii, constituent_d] = mAtlas->cmap().dartRanges().toLocalDart( c.dart() );
    const topology::Cell constituent_c( constituent_d, c.dim() );
    return mConstituents.at( patch_ii )->parentBasis( constituent_c );
}