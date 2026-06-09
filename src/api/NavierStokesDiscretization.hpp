#pragma once
#include <VectorConformingTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <VectorConformingHierarchicalTPSplineSpace.hpp>
#include <IndexOperations.hpp>

namespace api
{
    class NavierStokesDiscretization
    {
        public:

        virtual ~NavierStokesDiscretization() = default;

        virtual const Eigen::MatrixXd& geometryControlPoints() const = 0;
        virtual const Eigen::MatrixXd& controlPoints() const { return geometryControlPoints(); }
        virtual bool hasRationalGeometry() const { return false; }

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const = 0;

        virtual eval::SplineSpaceEvaluator& getH1() = 0;
        virtual eval::SplineSpaceEvaluator& getHDIV() = 0;
        virtual eval::SplineSpaceEvaluator& getL2() = 0;
        virtual eval::SplineSpaceEvaluator& getGeometry() { return getH1(); }

        virtual const eval::SplineSpaceEvaluator& getH1() const = 0;
        virtual const eval::SplineSpaceEvaluator& getHDIV() const = 0;
        virtual const eval::SplineSpaceEvaluator& getL2() const = 0;
        virtual const eval::SplineSpaceEvaluator& getGeometry() const { return getH1(); }
    };

    class NavierStokesTPDiscretization : public NavierStokesDiscretization
    {
        public:
        NavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                    const basis::KnotVector& kv_t,
                                    const size_t degree_s,
                                    const size_t degree_t,
                                    const Eigen::MatrixXd& cpts );

        virtual ~NavierStokesTPDiscretization() = default;

        virtual const Eigen::MatrixXd& geometryControlPoints() const override { return cpts; }

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const override { return cmap_bdry; }

        virtual eval::SplineSpaceEvaluator& getH1() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getHDIV() override { return HDIV; }
        virtual eval::SplineSpaceEvaluator& getL2() override { return L2; }

        virtual const eval::SplineSpaceEvaluator& getH1() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getHDIV() const override { return HDIV; }
        virtual const eval::SplineSpaceEvaluator& getL2() const override { return L2; }

        const basis::TPSplineSpace H1_ss;
        const basis::VectorConformingTPSplineSpace HDIV_ss;
        const basis::TPSplineSpace L2_ss;

        private:

        const topology::CombinatorialMapBoundary cmap_bdry;

        const Eigen::MatrixXd cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };

    class NURBSNavierStokesTPDiscretization : public NavierStokesTPDiscretization
    {
        public:
        NURBSNavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                           const basis::KnotVector& kv_t,
                                           const size_t degree_s,
                                           const size_t degree_t,
                                           const Eigen::MatrixXd& control_points,
                                           const Eigen::VectorXd& weights );

        NURBSNavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                           const basis::KnotVector& kv_t,
                                           const size_t degree_s,
                                           const size_t degree_t,
                                           const Eigen::MatrixXd& homogeneous_cpts );

        virtual ~NURBSNavierStokesTPDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const override { return euclidean_cpts; }
        virtual bool hasRationalGeometry() const override { return true; }

        virtual eval::SplineSpaceEvaluator& getGeometry() override { return Geometry; }
        virtual const eval::SplineSpaceEvaluator& getGeometry() const override { return Geometry; }

        const Eigen::VectorXd& weights() const { return cpt_weights; }

        private:
        const Eigen::MatrixXd euclidean_cpts;
        const Eigen::VectorXd cpt_weights;
        eval::NURBSSpaceEvaluator Geometry;
    };

    enum class PatchSide
    {
        S0,
        S1,
        T0,
        T1
    };

    class NavierStokesHierarchicalDiscretization : public NavierStokesDiscretization
    {
        public:
        NavierStokesHierarchicalDiscretization( const basis::KnotVector& kv_s,
                                                const basis::KnotVector& kv_t,
                                                const size_t degree_s,
                                                const size_t degree_t,
                                                const Eigen::MatrixXd& unrefined_cpts,
                                                const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine );

        virtual ~NavierStokesHierarchicalDiscretization() = default;

        virtual const Eigen::MatrixXd& geometryControlPoints() const override { return cpts; }

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const override { return cmap_bdry; }

        virtual eval::SplineSpaceEvaluator& getH1() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getHDIV() override { return HDIV; }
        virtual eval::SplineSpaceEvaluator& getL2() override { return L2; }

        virtual const eval::SplineSpaceEvaluator& getH1() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getHDIV() const override { return HDIV; }
        virtual const eval::SplineSpaceEvaluator& getL2() const override { return L2; }

        const basis::HierarchicalTPSplineSpace H1_ss;
        const basis::VectorConformingHierarchicalTPSplineSpace HDIV_ss;
        const basis::HierarchicalTPSplineSpace L2_ss;

        private:

        const topology::CombinatorialMapBoundary cmap_bdry;

        const Eigen::MatrixXd cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };
}
