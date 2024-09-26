#pragma once
#include <TPParametricAtlas.hpp>
#include <HierarchicalTPCombinatorialMap.hpp>

namespace param
{
    /// @brief A hierarchical parametric atlas over a series of nested tensor product parametric atlases.
    /// Tensor product atlases must be built over the refinement levels of the hierarchical topology.
    class HierarchicalTPParametricAtlas : public ParametricAtlas
    {
        public:
        HierarchicalTPParametricAtlas( const std::shared_ptr<const topology::HierarchicalTPCombinatorialMap>& cmap,
                                       const std::vector<std::shared_ptr<const TPParametricAtlas>>& refinement_levels );

        virtual ~HierarchicalTPParametricAtlas() = default;

        virtual const topology::HierarchicalTPCombinatorialMap& cmap() const override { return *mMap; }
        const std::shared_ptr<const topology::HierarchicalTPCombinatorialMap>& cmapPtr() const { return mMap; }

        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        const std::vector<std::shared_ptr<const TPParametricAtlas>>& refinementLevels() const { return mRefinementLevels; }

        private:
        const std::shared_ptr<const topology::HierarchicalTPCombinatorialMap> mMap;
        const std::vector<std::shared_ptr<const TPParametricAtlas>> mRefinementLevels;
    };
}