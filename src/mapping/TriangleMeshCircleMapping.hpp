#include <GeometricMapping.hpp>
#include <TriangleParametricAtlas.hpp>
#include <VertexPositionsFunc.hpp>
#include <map>

namespace mapping
{
    class TriangleMeshCircleMapping : public GeometricMapping
    {
        public:
        TriangleMeshCircleMapping( const param::TriangleParametricAtlas& atlas, const VertexPositionsFunc& vertex_positions );

        virtual const param::TriangleParametricAtlas& parametricAtlas() const override { return mAtlas; }

        virtual Eigen::VectorXd evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const override;

        private:
        const param::TriangleParametricAtlas& mAtlas;
        VertexPositionsFunc mPositions;
        std::map<size_t, double> mBoundaryAngles;
    };
}