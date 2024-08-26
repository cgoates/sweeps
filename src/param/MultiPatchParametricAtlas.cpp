#include <MultiPatchParametricAtlas.hpp>

namespace param
{
    MultiPatchParametricAtlas::MultiPatchParametricAtlas(
        const std::shared_ptr<const topology::MultiPatchCombinatorialMap>& cmap,
        const std::vector<std::shared_ptr<const TPParametricAtlas>>& constituents )
        : mCMap( cmap ), mSubAtlases( constituents )
    {}

    const topology::MultiPatchCombinatorialMap& MultiPatchParametricAtlas::cmap() const
    {
        return *mCMap;
    }

    const ParentDomain MultiPatchParametricAtlas::parentDomain( const topology::Cell& c ) const
    {
        const auto [patch_id, local_d] = mCMap->toLocalDart( c.dart() );
        const topology::Cell local_c( local_d, c.dim() );
        return mSubAtlases.at( patch_id )->parentDomain( local_c );
    }
    ParentPoint MultiPatchParametricAtlas::parentPoint( const topology::Vertex& v ) const
    {
        const auto [patch_id, local_d] = mCMap->toLocalDart( v.dart() );
        const topology::Vertex local_v( local_d );
        return mSubAtlases.at( patch_id )->parentPoint( local_v );
    }
    Vector6dMax MultiPatchParametricAtlas::parametricLengths( const topology::Cell& c ) const
    {
        const auto [patch_id, local_d] = mCMap->toLocalDart( c.dart() );
        const topology::Cell local_c( local_d, c.dim() );
        return mSubAtlases.at( patch_id )->parametricLengths( local_c );
    }
}