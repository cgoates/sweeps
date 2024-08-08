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
        TPBasisComplex( const std::shared_ptr<const param::TPParametricAtlas>& pa,
                        const std::shared_ptr<const BasisComplex>& source_complex,
                        const std::shared_ptr<const BasisComplex1d>& line_complex );
        virtual ~TPBasisComplex() = default;

        virtual const param::TPParametricAtlas& parametricAtlas() const override;
        const std::shared_ptr<const param::TPParametricAtlas>& parametricAtlasPtr() const { return mAtlas; }

        virtual ParentBasis parentBasis( const topology::Cell& ) const override;

        const BasisComplex& sourceComplex() const { return *mSourceComplex; }
        const BasisComplex1d& lineComplex() const { return *mLineComplex; }

        private:
        const std::shared_ptr<const param::TPParametricAtlas> mAtlas;
        const std::shared_ptr<const BasisComplex> mSourceComplex;
        const std::shared_ptr<const BasisComplex1d> mLineComplex;
    };
}