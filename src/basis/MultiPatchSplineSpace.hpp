#pragma once
#include <TPSplineSpace.hpp>
#include <MultiPatchBasisComplex.hpp>

namespace basis
{
    class MultiPatchSplineSpace : public SplineSpace
    {
        public:
        MultiPatchSplineSpace( const std::shared_ptr<const MultiPatchBasisComplex>& bc,
                               const std::vector<std::shared_ptr<const TPSplineSpace>>& constituents,
                               const std::vector<std::vector<FunctionId>>& func_ids );

        virtual ~MultiPatchSplineSpace() = default;

        virtual const MultiPatchBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        private:
        const std::shared_ptr<const MultiPatchBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const TPSplineSpace>> mSubSpaces;
        const std::vector<std::vector<FunctionId>> mFuncIds;
        size_t mNumFunctions;
    };
} // namespace basis