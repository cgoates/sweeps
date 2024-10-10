#pragma once
#include <Eigen/Dense>

namespace topology
{
    class Cell;
}

namespace param
{
    class ParametricAtlas;
    class ParentPoint;
}

namespace mapping
{
    class GeometricMapping
    {
        public:
        virtual const param::ParametricAtlas& parametricAtlas() const = 0;
        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const = 0;
        virtual size_t spatialDim() const = 0;
        /// NOTE: Eventually include a derivative call here
        /// NOTE: Perhaps use a localization call here eventually for performance reasons.
    };
}