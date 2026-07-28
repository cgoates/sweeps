#pragma once
#include <SplineSpace.hpp>
#include <ParentBasisEval.hpp>
#include <Cell.hpp>
#include <optional>
#include <CustomEigen.hpp>
#include <VertexPositionsFunc.hpp>
#include <vector>

namespace eval
{
    /// Evaluates a spline space on a localized element.
    ///
    /// Coordinate conventions:
    /// - Parent coordinates are element-local coordinates on the parent cell.
    /// - Parametric coordinates are coordinates in the spline's knot-space atlas.
    /// - Spatial coordinates are coordinates on the manifold defined by control points.
    ///
    /// Parent coordinates only locate an evaluation point inside the current element.
    /// Vector-conforming and top-form spline components remain expressed in the
    /// patch-parametric frame until an explicit spatial (Piola) transformation is
    /// applied.
    class SplineSpaceEvaluator
    {
        public:
        SplineSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv );

        virtual ~SplineSpaceEvaluator() = default;

        void localizeElement( const topology::Cell& );
        void localizeParentPoint( const param::ParentPoint& ppt );

        virtual Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& cpts ) const;
        /// Returns dx/ds with spatial components in rows and parent-coordinate
        /// directions in columns.
        virtual Eigen::MatrixXd evaluateParentToSpatialJacobian( const Eigen::MatrixXd& cpts ) const;
        /// Flattens the spatial Hessian into parent-coordinate columns ordered
        /// (ss, st, tt) in 2D and (ss, st, su, tt, tu, uu) in 3D.
        virtual Eigen::MatrixXd evaluateParentToSpatialHessian( const Eigen::MatrixXd& cpts ) const;

        /// Returns dx/dxi with spatial components in rows and parametric-coordinate
        /// directions in columns.
        virtual Eigen::MatrixXd evaluateParametricToSpatialJacobian( const Eigen::MatrixXd& cpts ) const;
        /// Flattens the spatial Hessian into parametric-coordinate columns ordered
        /// (xi-xi, xi-eta, eta-eta) in 2D and the analogous upper-triangular order in 3D.
        virtual Eigen::MatrixXd evaluateParametricToSpatialHessian( const Eigen::MatrixXd& cpts ) const;

        /// Evaluates the extracted spline components at the localized parent point.
        /// For vector-conforming and top-form spaces, the returned components are
        /// patch-parametric components.
        virtual Eigen::MatrixXd evaluateBasisValuesAtParentPoint() const;
        /// Columns are grouped first by parent-coordinate derivative and then by
        /// vector component.
        virtual Eigen::MatrixXd evaluateBasisFirstDerivativesWrtParentCoordinates() const;
        /// Columns use the symmetric parent-coordinate derivative ordering of
        /// evaluateParentToSpatialHessian, with vector components inside each group.
        virtual Eigen::MatrixXd evaluateBasisSecondDerivativesWrtParentCoordinates() const;

        /// Evaluates derivatives of the spline components with respect to
        /// patch-parametric coordinates. Columns are grouped first by parametric
        /// derivative direction and then by vector component.
        virtual Eigen::MatrixXd evaluateBasisFirstDerivativesWrtParametricCoordinates() const;

        virtual Eigen::VectorXd evaluateParametricPoint() const;

        size_t numDerivatives() const { return mNumDerivs; }

        const basis::SplineSpace& splineSpace() const { return mSpline; }

        protected:
        const size_t mNumDerivs;
        const basis::SplineSpace& mSpline;
        std::optional<topology::Cell> mCurrentCell;
        Vector6dMax mParametricLengths;
        Vector6dMax mParametricStart;
        std::vector<basis::FunctionId> mConnect;
        Eigen::MatrixXd mExOp;
        std::optional<ParentBasisEval> mLocalEvals;
        std::optional<Eigen::VectorXd> mLocalizedParentCoords;
    };

    /// Returns dJ_xi/dxi_k for every patch-parametric direction k.
    std::vector<Eigen::MatrixXd> evaluateParametricToSpatialJacobianFirstDerivatives(
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );

    /// Returns d(det J_xi)/dxi for a square parametric-to-spatial Jacobian.
    Eigen::VectorXd evaluateParametricToSpatialJacobianDeterminantFirstDerivatives(
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );

    /// Spatial H1 gradients, evaluated through the patch-parametric geometry map.
    Eigen::MatrixXd evaluateSpatialH1BasisFirstDerivatives(
        const SplineSpaceEvaluator& h1_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    /// Covariant Piola transformation from patch-parametric H(curl) components.
    Eigen::MatrixXd evaluateSpatialHCurlBasisValues(
        const SplineSpaceEvaluator& vec_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd evaluateSpatialHCurlBasisFirstDerivatives(
        const SplineSpaceEvaluator& vec_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    /// Contravariant Piola transformation from patch-parametric H(div) components.
    Eigen::MatrixXd evaluateSpatialHDivBasisValues(
        const SplineSpaceEvaluator& vec_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd evaluateSpatialHDivBasisFirstDerivatives(
        const SplineSpaceEvaluator& vec_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    /// Top-form Piola transformation from patch-parametric L2 coefficients.
    Eigen::MatrixXd evaluateSpatialL2BasisValues(
        const SplineSpaceEvaluator& l2_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );
    Eigen::MatrixXd evaluateSpatialL2BasisFirstDerivatives(
        const SplineSpaceEvaluator& l2_evals,
        const SplineSpaceEvaluator& geom_evals,
        const Eigen::MatrixXd& cpts );

    /// @brief Returns a VertexPositionsFunc by evaluating the manifold at each vertex.
    /// @param ss The spline space defining the manifold.
    /// @param cpts The control points defining the manifold, in a spatial_dim x n_control_pts sized matrix.
    VertexPositionsFunc vertexPositionsFromManifold( const basis::SplineSpace& ss, const Eigen::MatrixXd& cpts );


    class NURBSSpaceEvaluator : public SplineSpaceEvaluator
    {
        public:
        NURBSSpaceEvaluator( const basis::SplineSpace& ss, const size_t n_deriv )
            : SplineSpaceEvaluator( ss, n_deriv )
        {}

        virtual Eigen::VectorXd evaluateManifold( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateParentToSpatialJacobian( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateParentToSpatialHessian( const Eigen::MatrixXd& homogeneous_cpts ) const override;

        virtual Eigen::MatrixXd evaluateParametricToSpatialJacobian( const Eigen::MatrixXd& homogeneous_cpts ) const override;
        virtual Eigen::MatrixXd evaluateParametricToSpatialHessian( const Eigen::MatrixXd& homogeneous_cpts ) const override;

        virtual Eigen::MatrixXd evaluateBasisValuesAtParentPoint() const override;
        virtual Eigen::MatrixXd evaluateBasisFirstDerivativesWrtParentCoordinates() const override;
        virtual Eigen::MatrixXd evaluateBasisSecondDerivativesWrtParentCoordinates() const override;

        virtual Eigen::MatrixXd evaluateBasisFirstDerivativesWrtParametricCoordinates() const override;
    };
}
