#pragma once
#include <ParametricAtlas.hpp>
#include <CombinatorialMap1d.hpp>

namespace param
{
    /// @brief A parametric atlas over a 1d mesh.
    /// Defaults to unit parametric lengths.
    class ParametricAtlas1d : public ParametricAtlas
    {
        public:
        ParametricAtlas1d( const std::shared_ptr<const topology::CombinatorialMap1d>& cmap );
        ParametricAtlas1d( const std::shared_ptr<const topology::CombinatorialMap1d>& cmap,
                           const Eigen::VectorXd& lengths );
        virtual ~ParametricAtlas1d() = default;
        virtual const topology::CombinatorialMap1d& cmap() const override { return *mMap; }
        const std::shared_ptr<const topology::CombinatorialMap1d>& cmapPtr() const { return mMap; }
        virtual const ParentDomain parentDomain( const topology::Cell& c ) const override;
        virtual ParentPoint parentPoint( const topology::Vertex& v ) const override;
        virtual Vector6dMax parametricLengths( const topology::Cell& c ) const override;

        double totalLength() const { return mLengths.sum(); }

        private:
        const std::shared_ptr<const topology::CombinatorialMap1d> mMap;
        const ParentDomain mParentDomain = simplexDomain( 1 );
        Eigen::VectorXd mLengths;
    };

    std::vector<std::pair<topology::Cell, param::ParentPoint>> parentPointsOfParamPoints(
        const std::vector<double>& values, const param::ParametricAtlas1d& pa, const double param_tol );
}