#pragma once
#include <VectorConformingTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <VectorConformingHierarchicalTPSplineSpace.hpp>
#include <IndexOperations.hpp>

namespace assembly
{
    class DeRhamComplexDiscretization
    {
        public:

        virtual ~DeRhamComplexDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const = 0;

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const = 0;

        virtual eval::SplineSpaceEvaluator& getManifold() = 0;
        virtual eval::SplineSpaceEvaluator& getH1() = 0;
        virtual eval::SplineSpaceEvaluator& getHCURL() = 0;
        virtual eval::SplineSpaceEvaluator& getHDIV() = 0;
        virtual eval::SplineSpaceEvaluator& getL2() = 0;

        virtual const eval::SplineSpaceEvaluator& getManifold() const = 0;
        virtual const eval::SplineSpaceEvaluator& getH1() const = 0;
        virtual const eval::SplineSpaceEvaluator& getHCURL() const = 0;
        virtual const eval::SplineSpaceEvaluator& getHDIV() const = 0;
        virtual const eval::SplineSpaceEvaluator& getL2() const = 0;
    };

    class DeRhamComplexTPDiscretization : public DeRhamComplexDiscretization
    {
        public:
        DeRhamComplexTPDiscretization( const basis::KnotVector& kv_s,
                                    const basis::KnotVector& kv_t,
                                    const size_t degree_s,
                                    const size_t degree_t,
                                    const Eigen::Matrix2Xd& cpts );

        virtual ~DeRhamComplexTPDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const override { return cpts; }

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const override { return cmap_bdry; }

        virtual eval::SplineSpaceEvaluator& getManifold() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getH1() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getHCURL() override { return HCURL; }
        virtual eval::SplineSpaceEvaluator& getHDIV() override { return HDIV; }
        virtual eval::SplineSpaceEvaluator& getL2() override { return L2; }

        virtual const eval::SplineSpaceEvaluator& getManifold() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getH1() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getHCURL() const override { return HCURL; }
        virtual const eval::SplineSpaceEvaluator& getHDIV() const override { return HDIV; }
        virtual const eval::SplineSpaceEvaluator& getL2() const override { return L2; }

        // TODO: Add optional Manifold spline
        const basis::TPSplineSpace H1_ss;
        const basis::VectorConformingTPSplineSpace HCURL_ss;
        const basis::VectorConformingTPSplineSpace HDIV_ss;
        const basis::TPSplineSpace L2_ss;

        private:

        const topology::CombinatorialMapBoundary cmap_bdry;

        const Eigen::MatrixXd cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HCURL;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };

    class DeRhamComplexHierarchicalDiscretization : public DeRhamComplexDiscretization
    {
        public:
        DeRhamComplexHierarchicalDiscretization( const basis::KnotVector& kv_s,
                                                const basis::KnotVector& kv_t,
                                                const size_t degree_s,
                                                const size_t degree_t,
                                                const Eigen::Matrix2Xd& unrefined_cpts,
                                                const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine );

        DeRhamComplexHierarchicalDiscretization(
            const basis::KnotVector& kv_s,
            const basis::KnotVector& kv_t,
            const basis::KnotVector& kv_u,
            const size_t degree_s,
            const size_t degree_t,
            const size_t degree_u,
            const Eigen::Matrix3Xd& unrefined_cpts,
            const std::vector<std::vector<std::array<size_t, 3>>>& elems_to_refine );

        virtual ~DeRhamComplexHierarchicalDiscretization() = default;

        virtual const Eigen::MatrixXd& controlPoints() const override { return cpts; }

        virtual const topology::CombinatorialMapBoundary& cmapBdry() const override { return cmap_bdry; }

        virtual eval::SplineSpaceEvaluator& getManifold() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getH1() override { return H1; }
        virtual eval::SplineSpaceEvaluator& getHCURL() override { return HCURL; }
        virtual eval::SplineSpaceEvaluator& getHDIV() override { return HDIV; }
        virtual eval::SplineSpaceEvaluator& getL2() override { return L2; }

        virtual const eval::SplineSpaceEvaluator& getManifold() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getH1() const override { return H1; }
        virtual const eval::SplineSpaceEvaluator& getHCURL() const override { return HCURL; }
        virtual const eval::SplineSpaceEvaluator& getHDIV() const override { return HDIV; }
        virtual const eval::SplineSpaceEvaluator& getL2() const override { return L2; }

        // TODO: Add optional Manifold spline
        const basis::HierarchicalTPSplineSpace H1_ss;
        const basis::VectorConformingHierarchicalTPSplineSpace HCURL_ss;
        const basis::VectorConformingHierarchicalTPSplineSpace HDIV_ss;
        const basis::HierarchicalTPSplineSpace L2_ss;

        private:

        const topology::CombinatorialMapBoundary cmap_bdry;

        const Eigen::MatrixXd cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HCURL;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };

    void localizeElement( DeRhamComplexDiscretization& nsd, const topology::Cell& elem );

    void localizePoint( DeRhamComplexDiscretization& nsd, const param::ParentPoint& ppt );
}