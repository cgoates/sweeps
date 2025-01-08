#pragma once
#include <Eigen/Dense>
#include <optional>
#include <Cell.hpp>
#include <ParentPoint.hpp>

namespace param
{
    class ParametricAtlas;
}

namespace mapping
{
    class GeometricMapping
    {
        public:
        virtual const param::ParametricAtlas& parametricAtlas() const = 0;
        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const = 0;
        virtual size_t spatialDim() const = 0;
        virtual std::optional<std::pair<topology::Cell, param::ParentPoint>> maybeInverse( const Eigen::Vector2d& pt ) const = 0;
        virtual std::pair<topology::Cell, param::ParentPoint> closestPoint( const Eigen::Vector3d& ) const
        {
            throw std::runtime_error( "closestPoint not implemented for this mapping" );
        }
        /// NOTE: Eventually include a derivative call here
        /// NOTE: Perhaps use a localization call here eventually for performance reasons.
    };
}