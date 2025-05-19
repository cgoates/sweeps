#pragma once
#include <GeometricMapping.hpp>
#include <TriangleParametricAtlas.hpp>
#include <VertexPositionsFunc.hpp>
#include <map>

namespace mapping
{
    class AABB
    {
        public:
        AABB( const Vector3dMax& min, const Vector3dMax& max ) : mMin( min ), mMax( max ) {}

        const Vector3dMax& min() const { return mMin; }
        const Vector3dMax& max() const { return mMax; }

        bool contains( const Vector3dMax& pt ) const
        {
            return ( pt.array() >= mMin.array() ).all() and ( pt.array() <= mMax.array() ).all();
        }

        private:
        Vector3dMax mMin;
        Vector3dMax mMax;
    };

    std::ostream& operator<<( std::ostream& os, const AABB& aabb );

    class TriangleMeshMapping : public GeometricMapping
    {
        public:
        TriangleMeshMapping( const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
                             const VertexPositionsFunc& vertex_positions,
                             const size_t dim,
                             const bool build_bounding_boxes = true );

        virtual ~TriangleMeshMapping() = default;

        virtual const param::TriangleParametricAtlas& parametricAtlas() const override { return *mAtlas; }

        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const override;

        virtual size_t spatialDim() const override { return mDim; }

        std::optional<param::ParentPoint> maybeInverse( const topology::Face& f, const Vector3dMax& pt ) const;

        //FIXME: PROBABLY VERY SLOW
        virtual std::optional<std::pair<topology::Cell, param::ParentPoint>> maybeInverse( const Vector3dMax& pt ) const override;

        virtual std::pair<topology::Cell, param::ParentPoint> closestPoint( const Vector3dMax& pt ) const override;

        virtual const VertexPositionsFunc& vertPositions() const override { return mPositions; }

        private:
        size_t vertexIndex( const topology::Vertex& v ) const;
        const std::shared_ptr<const param::TriangleParametricAtlas> mAtlas;
        VertexPositionsFunc mPositions;
        std::map<topology::Face, AABB> mBoundingBoxes;
        size_t mDim;
    };
}