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
        DivConfBasisComplex( const std::shared_ptr<const BasisComplex>& primal_complex );
        virtual ~DivConfBasisComplex() = default;

        virtual const param::ParametricAtlas& parametricAtlas() const override;

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        private:
        const std::shared_ptr<const BasisComplex> mPrimalComplex;
    };
}