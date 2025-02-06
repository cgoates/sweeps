#include <HierarchicalMultiPatchParametricAtlas.hpp>

using namespace param;

std::vector<std::shared_ptr<const HierarchicalTPParametricAtlas>>
    initializeConstituents( const topology::HierarchicalMultiPatchCombinatorialMap& cmap,
                            const std::vector<std::shared_ptr<const MultiPatchParametricAtlas>>& refinement_levels )
{
    if( cmap.numLevels() != refinement_levels.size() )
        throw std::invalid_argument( "Inconsistent refinement levels and cmap provided to HierarchicalMultiPatchParametricAtlas." );

    std::vector<std::vector<std::shared_ptr<const TPParametricAtlas>>> constituent_levels( cmap.constituents().size() );
    for( const auto& level : refinement_levels )
    {
        if( level->constituents().size() != cmap.constituents().size() )
            throw std::invalid_argument( "All refinement levels must have the same number of constituents." );
        for( size_t patch_ii = 0; patch_ii < cmap.constituents().size(); patch_ii++ )
        {
            constituent_levels.at( patch_ii ).push_back( level->constituents().at( patch_ii ) );
        }
    }

    std::vector<std::shared_ptr<const HierarchicalTPParametricAtlas>> out;
    out.reserve( cmap.constituents().size() );
    for( size_t i = 0; i < cmap.constituents().size(); i++ )
    {
        out.push_back( std::make_shared<HierarchicalTPParametricAtlas>( cmap.constituents().at( i ), constituent_levels.at( i ) ) );
    }
    return out;
}

HierarchicalMultiPatchParametricAtlas::HierarchicalMultiPatchParametricAtlas( const std::shared_ptr<const topology::HierarchicalMultiPatchCombinatorialMap>& cmap,
                                    const std::vector<std::shared_ptr<const MultiPatchParametricAtlas>>& refinement_levels )
    : mMap( cmap ), mRefinementLevels( refinement_levels ), mConstituents( initializeConstituents( *cmap, refinement_levels ) )
{}

const ParentDomain HierarchicalMultiPatchParametricAtlas::parentDomain( const topology::Cell& c ) const
{
    const auto [patch_ii, constituent_d] = mMap->dartRanges().toLocalDart( c.dart() );
    const topology::Cell constituent_c( constituent_d, c.dim() );
    return mConstituents.at( patch_ii )->parentDomain( constituent_c );
}

ParentPoint HierarchicalMultiPatchParametricAtlas::parentPoint( const topology::Vertex& v ) const
{
    const auto [patch_ii, constituent_d] = mMap->dartRanges().toLocalDart( v.dart() );
    const topology::Vertex constituent_v( constituent_d );
    return mConstituents.at( patch_ii )->parentPoint( constituent_v );
}
Vector6dMax HierarchicalMultiPatchParametricAtlas::parametricLengths( const topology::Cell& c ) const
{
    const auto [patch_ii, constituent_d] = mMap->dartRanges().toLocalDart( c.dart() );
    const topology::Cell constituent_c( constituent_d, c.dim() );
    return mConstituents.at( patch_ii )->parametricLengths( constituent_c );
}