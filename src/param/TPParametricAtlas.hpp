#pragma once
#include <ParametricAtlas.hpp>
#include <ParametricAtlas1d.hpp>
#include <TPCombinatorialMap.hpp>

namespace param
{
    /// @brief A parametric atlas over a tensor product mesh. Source and line atlases must
    /// be built over the cmaps that are source and line cmaps of the tensor product cmap.
    class TPParametricAtlas : public ParametricAtlas
    {
        public:
        TPParametricAtlas( const std::shared_ptr<const topology::TPCombinatorialMap>& cmap,
                           const std::shared_ptr<const ParametricAtlas>& source_atlas,
                           const std::shared_ptr<const ParametricAtlas1d>& line_atlas );
        virtual ~TPParametricAtlas() = default;
        virtual const topology::TPCombinatorialMap& cmap() const override { return *mMap; }
        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        const ParametricAtlas& source() const { return *mSourceParam; }
        const ParametricAtlas1d& line() const { return *mLineParam; }

        private:
        const std::shared_ptr<const topology::TPCombinatorialMap> mMap;
        const std::shared_ptr<const ParametricAtlas> mSourceParam;
        const std::shared_ptr<const ParametricAtlas1d> mLineParam;
    };
}