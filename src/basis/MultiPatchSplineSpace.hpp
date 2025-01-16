#pragma once
#include <TPSplineSpace.hpp>
#include <MultiPatchBasisComplex.hpp>

namespace basis
{
    class MultiPatchSplineSpace : public SplineSpace
    {
        public:
        MultiPatchSplineSpace( const std::shared_ptr<const MultiPatchBasisComplex>& bc,
                               const std::vector<std::shared_ptr<const TPSplineSpace>>& constituents );

        virtual ~MultiPatchSplineSpace() = default;

        virtual const MultiPatchBasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const std::vector<std::vector<FunctionId>>& functionIdMap() const { return mFuncIds; }

        const std::vector<std::shared_ptr<const TPSplineSpace>>& subSpaces() const { return mSubSpaces; }

        private:
        const std::shared_ptr<const MultiPatchBasisComplex> mBasisComplex;
        const std::vector<std::shared_ptr<const TPSplineSpace>> mSubSpaces;
        std::vector<std::vector<FunctionId>> mFuncIds;
        size_t mNumFunctions;
    };

    MultiPatchSplineSpace buildMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const std::map<std::pair<size_t, topology::Dart>, std::pair<size_t, topology::Dart>>& connections );

    Eigen::MatrixX3d multiPatchCoefficients( const MultiPatchSplineSpace& ss,
                                             const std::vector<Eigen::MatrixX3d>& patch_coeffs );
} // namespace basis