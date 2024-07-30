#pragma once
#include <SplineSpace.hpp>
#include <ParentBasisEval.hpp>
#include <Cell.hpp>

namespace eval
{
    class SplineSpaceEvaluator
    {
        public:
        SplineSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv );

        void localizeElement( const topology::Cell& );
        void localizePoint( const param::ParentPoint& ppt );

        Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& cpts ) const;
        Eigen::MatrixXd evaluateJacobian( const Eigen::MatrixXd& cpts ) const;
        Eigen::MatrixXd evaluateHessian( const Eigen::MatrixXd& cpts ) const;// FIXME: Document the returned matrix ordering. Flattening a 3-tensor into a matrix here.
        Eigen::MatrixXd evaluatePiola( const Eigen::MatrixXd& cpts ) const;

        Eigen::MatrixXd evaluateBasis() const;
        Eigen::MatrixXd evaluateFirstDerivatives() const;
        Eigen::MatrixXd evaluateSecondDerivatives() const;// FIXME: Document the returned matrix ordering.

        size_t numDerivatives() const { return mNumDerivs; }

        private:
        const size_t mNumDerivs;
        const basis::SplineSpace& mSpline;
        topology::Cell mCurrentCell;
        std::vector<basis::FunctionId> mConnect;
        Eigen::MatrixXd mExOp;
        ParentBasisEval mLocalEvals;
    };

    Eigen::MatrixXd piolaTransformedVectorBasis( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedVectorFirstDerivatives( const SplineSpaceEvaluator& vec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::VectorXd piolaTransformedBivectorBasis( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd piolaTransformedBivectorFirstDerivatives( const SplineSpaceEvaluator& bivec_evals, const SplineSpaceEvaluator& geom_evals, const Eigen::MatrixXd& cpts );
}