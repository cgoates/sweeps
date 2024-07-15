#pragma once
#include <SplineSpace.hpp>
#include <BasisComplex1d.hpp>
#include <KnotVector.hpp>
#include <map>

namespace basis
{
    class BasisComplex1d;
    class KnotVector;

    /// @brief A spline space that is encoded simply by storing the extraction operators and connectivity
    /// for each element in the BasisComplex.
    class BSplineSpace1d : public SplineSpace
    {
        public:
        BSplineSpace1d( const BasisComplex1d& bc, const KnotVector& kv );
        virtual ~BSplineSpace1d() = default;

        virtual const BasisComplex1d& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        virtual size_t numFunctions() const override;

        const KnotVector& knotVector() const { return mKnotVector; }

        private:
        const BasisComplex1d& mBasisComplex;
        const KnotVector mKnotVector;
        std::map<size_t, Eigen::MatrixXd> mExtractionOps;
        std::map<size_t, std::vector<FunctionId>> mConnectivity;
    };
}