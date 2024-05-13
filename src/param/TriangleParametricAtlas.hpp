#pragma once
#include <ParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>

namespace param
{
    /// @brief A parametric atlas over a triangle mesh.  Providing a combinatorial map
    /// to this class that is not triangle-only will result in errors.
    class TriangleParametricAtlas : public ParametricAtlas
    {
        public:
        TriangleParametricAtlas( const topology::CombinatorialMap& cmap ) : mMap( cmap ) {}
        virtual const topology::CombinatorialMap& cmap() const override { return mMap; }
        virtual const ParentDomain& parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;

        private:
        const topology::CombinatorialMap& mMap;
        const ParentDomain mParentDomain = simplexDomain( 2 );
        const std::array<Eigen::Vector2d, 3> mPoints = { Eigen::Vector2d{1, 1}, Eigen::Vector2d{0, 1}, Eigen::Vector2d{1, 0} };
    };
}