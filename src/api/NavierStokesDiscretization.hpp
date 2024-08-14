#pragma once
#include <DivConfTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <CombinatorialMapBoundary.hpp>

namespace api
{
    struct NavierStokesDiscretization
    {
        NavierStokesDiscretization( const basis::KnotVector& kv_s,
                                    const basis::KnotVector& kv_t,
                                    const size_t degree_s,
                                    const size_t degree_t,
                                    const Eigen::MatrixX2d& cpts );

        const basis::TPSplineSpace H1_ss;
        const basis::DivConfTPSplineSpace HDIV_ss;
        const basis::TPSplineSpace L2_ss;

        const topology::CombinatorialMapBoundary cmap_bdry;

        const Eigen::MatrixX2d cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };

    enum class PatchSide
    {
        S0,
        S1,
        T0,
        T1
    };
}