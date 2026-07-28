#pragma once
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>
#include <ParentBasis.hpp>

namespace basis
{
    enum class ConformingType
    {
        Divergence,
        Curl
    };
    /// @brief A parent basis complex that provides div or curl conforming parent bases.
    /// Layers over any basis complex that only includes cube-like cells.
    class VectorConformingBasisComplex : public BasisComplex
    {
        public:
        VectorConformingBasisComplex( const std::shared_ptr<const BasisComplex>& primal_complex,
                                      const ConformingType conforming_type = ConformingType::Divergence );
        virtual ~VectorConformingBasisComplex() = default;

        virtual const param::ParametricAtlas& parametricAtlas() const override;

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        ConformingType conformingType() const { return mConformingType; }
        const std::shared_ptr<const BasisComplex>& primalComplexPtr() const { return mPrimalComplex; }

        private:
        const std::shared_ptr<const BasisComplex> mPrimalComplex;
        const ConformingType mConformingType;
    };
}
