#pragma once
#include <MultiPatchParametricAtlas.hpp>
#include <HierarchicalMultiPatchCombinatorialMap.hpp>
#include <HierarchicalTPParametricAtlas.hpp>

namespace param
{
    /// @brief A hierarchical parametric atlas over a series of nested tensor product parametric atlases.
    /// Tensor product atlases must be built over the refinement levels of the hierarchical topology.
    class HierarchicalMultiPatchParametricAtlas : public ParametricAtlas
    {
        public:
        HierarchicalMultiPatchParametricAtlas(
            const std::shared_ptr<const topology::HierarchicalMultiPatchCombinatorialMap>& cmap,
            const std::vector<std::shared_ptr<const MultiPatchParametricAtlas>>& refinement_levels );

        virtual ~HierarchicalMultiPatchParametricAtlas() = default;

        virtual const topology::HierarchicalMultiPatchCombinatorialMap& cmap() const override { return *mMap; }
        const std::shared_ptr<const topology::HierarchicalMultiPatchCombinatorialMap>& cmapPtr() const { return mMap; }

        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        const std::vector<std::shared_ptr<const MultiPatchParametricAtlas>>& refinementLevels() const { return mRefinementLevels; }

        private:
        const std::shared_ptr<const topology::HierarchicalMultiPatchCombinatorialMap> mMap;
        const std::vector<std::shared_ptr<const MultiPatchParametricAtlas>> mRefinementLevels;
        const std::vector<std::shared_ptr<const HierarchicalTPParametricAtlas>> mConstituents;
    };
}