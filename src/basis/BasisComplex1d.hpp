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
        BasisComplex1d( const param::ParametricAtlas1d& pa, const uint degree );
        ~BasisComplex1d() = default;

        virtual const param::ParametricAtlas1d& parametricAtlas() const override;

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const param::ParametricAtlas1d& mAtlas;
        const ParentBasis mParentBasis;
    };
}