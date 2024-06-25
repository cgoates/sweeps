#pragma once
#include <BasisComplex.hpp>
#include <TPParametricAtlas.hpp>
#include <ParentBasis.hpp>
#include <BasisComplex1d.hpp>

namespace basis
{
    /// @brief A parent basis complex over a ParametricAtlas1d.
    /// Supports only uniform degree for now.
    class TPBasisComplex : public BasisComplex
    {
        public:
        TPBasisComplex( const param::TPParametricAtlas& pa, const BasisComplex& source_complex, const BasisComplex1d& line_complex );
        ~TPBasisComplex() = default;

        virtual const param::TPParametricAtlas& parametricAtlas() const override;

        virtual const ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const param::TPParametricAtlas& mAtlas;
        const BasisComplex& mSourceComplex;
        const BasisComplex1d& mLineComplex;
    };
}