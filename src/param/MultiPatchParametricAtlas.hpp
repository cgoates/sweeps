#pragma once
#include <TPParametricAtlas.hpp>
#include <MultiPatchCombinatorialMap.hpp>

namespace param
{
    class MultiPatchParametricAtlas : public ParametricAtlas
    {
        public:
        MultiPatchParametricAtlas( const std::shared_ptr<const topology::MultiPatchCombinatorialMap>& cmap,
                                   const std::vector<std::shared_ptr<const TPParametricAtlas>>& constituents );

        virtual ~MultiPatchParametricAtlas() = default;

        virtual const topology::MultiPatchCombinatorialMap& cmap() const override;
        const std::shared_ptr<const topology::MultiPatchCombinatorialMap>& cmapPtr() const { return mCMap; }
        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        const std::vector<std::shared_ptr<const TPParametricAtlas>>& constituents() const { return mSubAtlases; }

        private:
        const std::shared_ptr<const topology::MultiPatchCombinatorialMap> mCMap;
        const std::vector<std::shared_ptr<const TPParametricAtlas>> mSubAtlases;
    };
}