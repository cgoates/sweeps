#pragma once
#include <ParametricAtlas.hpp>

namespace param
{
    /// @brief A parametric atlas over a triangle mesh.  Providing a combinatorial map
    /// to this class that is not triangle-only will result in errors.
    class TriangleParametricAtlas : public ParametricAtlas
    {
        public:
        TriangleParametricAtlas( const std::shared_ptr<const topology::CombinatorialMap>& cmap ) : mMap( cmap ) {}
        virtual ~TriangleParametricAtlas() = default;
        virtual const topology::CombinatorialMap& cmap() const override { return *mMap; }
        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        private:
        const std::shared_ptr<const topology::CombinatorialMap> mMap;
        const ParentDomain mParentDomain = simplexDomain( 2 );
    };
}