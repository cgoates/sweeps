#pragma once
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>
#include <ParentBasis.hpp>

namespace basis
{
    /// @brief A parent basis complex that provides div conforming parent bases.
    /// Layers over any basis complex that only includes cube-like cells.
    class DivConfBasisComplex : public BasisComplex
    {
        public:
        DivConfBasisComplex( const BasisComplex& primal_complex );
        ~DivConfBasisComplex() = default;

        virtual const param::ParametricAtlas& parametricAtlas() const override;

        virtual const ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const BasisComplex& mPrimalComplex;
    };
}