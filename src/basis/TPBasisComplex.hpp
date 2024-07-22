#pragma once
#include <BasisComplex.hpp>
#include <TPParametricAtlas.hpp>
#include <ParentBasis.hpp>
#include <BasisComplex1d.hpp>

namespace basis
{
    /// @brief A parent basis complex over a TPParametricAtlas.
    ///
    class TPBasisComplex : public BasisComplex
    {
        public:
        TPBasisComplex( const param::TPParametricAtlas& pa, const BasisComplex& source_complex, const BasisComplex1d& line_complex );
        virtual ~TPBasisComplex() = default;

        virtual const param::TPParametricAtlas& parametricAtlas() const override;

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        const BasisComplex& sourceComplex() const { return mSourceComplex; }
        const BasisComplex1d& lineComplex() const { return mLineComplex; }

        private:
        const param::TPParametricAtlas& mAtlas;
        const BasisComplex& mSourceComplex;
        const BasisComplex1d& mLineComplex;
    };
}