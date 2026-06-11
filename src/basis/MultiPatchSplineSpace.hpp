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
                               const bool connect_interfaces = true );

        virtual ~MultiPatchSplineSpace() = default;

        virtual const MultiPatchBasisComplex& basisComplex() const override;
        const std::shared_ptr<const MultiPatchBasisComplex>& basisComplexPtr() const { return mBasisComplex; }

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

    MultiPatchSplineSpace buildMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const topology::MultiPatchCombinatorialMap::InternalConnectionsMap& connections );

    MultiPatchSplineSpace buildDiscontinuousMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const std::map<std::pair<size_t, topology::Dart>, std::pair<size_t, topology::Dart>>& connections );

    MultiPatchSplineSpace buildDiscontinuousMultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const topology::MultiPatchCombinatorialMap::InternalConnectionsMap& connections );

    struct DegreeAndKnotVector
    {
        SmallVector<size_t, 3> degrees;
        SmallVector<basis::KnotVector, 3> kvs;
    };

    MultiPatchSplineSpace degreeRefineOrCoarsen( const MultiPatchSplineSpace& ss,
                                                 const std::function<DegreeAndKnotVector( const size_t )>& degree_and_kv_func );

    Eigen::MatrixXd multiPatchCoefficients( const MultiPatchSplineSpace& ss,
                                            const std::vector<Eigen::MatrixXd>& patch_coeffs );

    std::vector<Eigen::MatrixXd> splitMultiPatchCoefficients( const MultiPatchSplineSpace& ss,
                                                              const Eigen::MatrixXd& global_coeffs );

    /// A helper function to get iteration variables for traversing the control points of a constituent TPSplineSpace
    /// This allows for connection of control points between the patches.
    std::tuple<util::IndexVec, util::IndexVec, SmallVector<std::variant<bool, size_t>, 3>>
        getIterVars( const TPSplineSpace& constituent, const topology::Cell& corner, const bool reverse_dart = false );
} // namespace basis
