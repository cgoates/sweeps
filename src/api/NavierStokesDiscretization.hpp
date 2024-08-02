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

        const topology::CombinatorialMap1d cmap_1d_s;
        const topology::CombinatorialMap1d cmap_1d_t;
        const param::ParametricAtlas1d param_1d_s;
        const param::ParametricAtlas1d param_1d_t;
        const basis::BasisComplex1d bc_1d_s;
        const basis::BasisComplex1d bc_1d_t;
        const basis::BSplineSpace1d ss_1d_s;
        const basis::BSplineSpace1d ss_1d_t;

        const topology::TPCombinatorialMap cmap_2d;
        const topology::CombinatorialMapBoundary cmap_bdry;
        const param::TPParametricAtlas param_2d;
        const basis::TPBasisComplex H1_bc;
        const basis::TPSplineSpace H1_ss;

        const basis::DivConfBasisComplex HDIV_bc;
        const basis::DivConfTPSplineSpace HDIV_ss;

        const basis::TPBasisComplex L2_bc;
        const basis::TPSplineSpace L2_ss;

        const Eigen::MatrixX2d cpts;

        eval::SplineSpaceEvaluator H1;
        eval::SplineSpaceEvaluator HDIV;
        eval::SplineSpaceEvaluator L2;
    };
}