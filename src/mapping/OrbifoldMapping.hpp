#pragma once
#include <GeometricMapping.hpp>
#include <TriangleMeshMapping.hpp>
#include <VertexPositionsFunc.hpp>
#include <SimplexUtilities.hpp>
#include <AABB.hpp>
#include <map>

namespace mapping
{
    struct OrbifoldTriangle
    {
        topology::Face face;
        SmallVector<std::pair<AABB, Triangle<2>>, 2> triangles;
    };

    class OrbifoldMapping : public GeometricMapping
    {
        public:
        OrbifoldMapping( const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
                         const VertexPositionsFunc& vertex_positions );

        virtual ~OrbifoldMapping() = default;

        virtual const param::TriangleParametricAtlas& parametricAtlas() const override;

        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const override;

        virtual size_t spatialDim() const override;

        /*
         * From the figure of the David bust orbifold, we can see that transforming the point by an integer distance in
         * x or y rotates it by +-90 degrees if the integer is odd, and does not rotate it if the integer is even. So to
         * get the point into the canonical domain, just take the non-integer part for its coordinates, and then rotate
         * as needed.  Note that this assumes the canonical domain is [0,1]x[0,1], which is not how I originally coded
         * the orbifold embedding. FIX THIS
         */
        virtual std::optional<std::pair<topology::Cell, param::ParentPoint>> maybeInverse( const Vector3dMax& pt ) const override;

        virtual std::pair<topology::Cell, param::ParentPoint> closestPoint( const Vector3dMax& pt ) const override;

        virtual const VertexPositionsFunc& vertPositions() const override;

        void iterateTriangles( const std::function<void( const Triangle<2>& )>& callback ) const
        {
            for( const auto& tri : mTriangles )
            {
                for( const auto& [bb, t] : tri.triangles )
                {
                    callback( t );
                }
            }
        }

        private:
        const TriangleMeshMapping mTriangleMapping;
        std::vector<OrbifoldTriangle> mTriangles;
    };
}