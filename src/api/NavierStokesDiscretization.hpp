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

        virtual const Eigen::MatrixXd& controlPoints() const = 0;

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const = 0;

        virtual eval::SplineSpaceEvaluator& getH1() = 0;
        virtual eval::SplineSpaceEvaluator& getHDIV() = 0;
        virtual eval::SplineSpaceEvaluator& getL2() = 0;

        virtual const eval::SplineSpaceEvaluator& getH1() const = 0;
        virtual const eval::SplineSpaceEvaluator& getHDIV() const = 0;
        virtual const eval::SplineSpaceEvaluator& getL2() const = 0;
    };

    class NavierStokesTPDiscretization : public NavierStokesDiscretization
    {
        public:
        NavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                    const basis::KnotVector& kv_t,
                                    const size_t degree_s,
                                    const size_t degree_t,
                                    const Eigen::Matrix2Xd& cpts );

        virtual ~NavierStokesTPDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const override { return cpts; }

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

    enum class PatchSide
    {
        S0,
        S1,
        T0,
        T1,
        U0,
        U1
    };

    class NavierStokesHierarchicalDiscretization : public NavierStokesDiscretization
    {
        public:
        NavierStokesHierarchicalDiscretization( const basis::KnotVector& kv_s,
                                                const basis::KnotVector& kv_t,
                                                const size_t degree_s,
                                                const size_t degree_t,
                                                const Eigen::Matrix2Xd& unrefined_cpts,
                                                const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine );

        virtual ~NavierStokesHierarchicalDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const override { return cpts; }

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

    void localizeElement( api::NavierStokesDiscretization& nsd, const topology::Cell& elem );

    void localizePoint( api::NavierStokesDiscretization& nsd, const param::ParentPoint& ppt );
}