#pragma once
#include <SplineSpace.hpp>
#include <BasisComplex.hpp>
#include <map>

namespace basis
{
    class BasisComplex1d;
    class KnotVector;

    /// @brief A spline space that is encoded simply by storing the extraction operators and connectivity
    /// for each element in the BasisComplex.
    class GenericSplineSpace : public SplineSpace
    {
        public:
        GenericSplineSpace( const BasisComplex& bc,
                            const std::map<size_t, Eigen::MatrixXd>& ex_ops,
                            const std::map<size_t, std::vector<FunctionId>>& connectivity );

        virtual const BasisComplex& basisComplex() const override;

        virtual Eigen::MatrixXd extractionOperator( const topology::Cell& ) const override;

        virtual std::vector<FunctionId> connectivity( const topology::Cell& ) const override;

        private:
        const BasisComplex& mBasisComplex;
        std::map<size_t, Eigen::MatrixXd> mExtractionOps;
        std::map<size_t, std::vector<FunctionId>> mConnectivity;
    };

    GenericSplineSpace coxDeBoor( const BasisComplex1d& bc, const KnotVector& kv );
}