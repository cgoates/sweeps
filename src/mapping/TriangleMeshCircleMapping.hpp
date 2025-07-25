#include <GeometricMapping.hpp>
#include <TriangleParametricAtlas.hpp>
#include <VertexPositionsFunc.hpp>
#include <TriangleMeshMapping.hpp>
#include <map>

namespace mapping
{
    class TriangleMeshCircleMapping : public GeometricMapping
    {
        public:
        TriangleMeshCircleMapping( const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
                                   const VertexPositionsFunc& vertex_positions );

        virtual ~TriangleMeshCircleMapping() = default;

        virtual const param::TriangleParametricAtlas& parametricAtlas() const override { return *mAtlas; }

        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const override;

        virtual size_t spatialDim() const override { return 2; }

        std::optional<param::ParentPoint> maybeInverse( const topology::Face& f, const Eigen::Vector2d& pt ) const;

        virtual const VertexPositionsFunc& vertPositions() const override { return mTriMapping.vertPositions(); }

        //FIXME: PROBABLY VERY SLOW
        virtual std::optional<std::pair<topology::Cell, param::ParentPoint>> maybeInverse( const Vector3dMax& pt ) const override;

        private:
        const std::shared_ptr<const param::TriangleParametricAtlas> mAtlas;
        const TriangleMeshMapping mTriMapping;
        std::map<size_t, double> mBoundaryAngles;
        std::map<topology::Face, AABB> mBoundingBoxes;
    };
}