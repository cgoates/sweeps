#pragma once
#include <BasisComplex.hpp>
#include <ParametricAtlas1d.hpp>
#include <ParentBasis.hpp>

namespace basis
{
    /// @brief A parent basis complex over a ParametricAtlas1d.
    /// Supports only uniform degree for now.
    class BasisComplex1d : public BasisComplex
    {
        public:
        BasisComplex1d( const std::shared_ptr<const param::ParametricAtlas1d>& pa, const uint degree );
        virtual ~BasisComplex1d() = default;

        virtual const param::ParametricAtlas1d& parametricAtlas() const override;
        const std::shared_ptr<const param::ParametricAtlas1d>& parametricAtlasPtr() const { return mAtlas; }

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        ParentBasis defaultParentBasis() const { return mParentBasis; }

        private:
        const std::shared_ptr<const param::ParametricAtlas1d> mAtlas;
        const ParentBasis mParentBasis;
    };

    BasisComplex1d reduceDegree( const BasisComplex1d& bc );
}