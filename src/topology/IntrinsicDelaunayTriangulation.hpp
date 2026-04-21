#pragma once
#include <CombinatorialMap.hpp>
#include <unordered_map>
#include <vector>
#include <VertexPositionsFunc.hpp>
#include <memory>
#include <Eigen/Core>

namespace topology
{
    /// @brief A combinatorial map that represents the intrinsic Delaunay triangulation of a
    /// triangle mesh embedded in 3D. Starting from a given triangle mesh, interior edges are
    /// iteratively flipped until every interior edge satisfies the Delaunay condition
    /// (sum of opposite angles ≤ π). phi2 is unchanged; only phi1/phi(-1) are overridden.
    class IntrinsicDelaunayTriangulation : public CombinatorialMap
    {
        public:
        IntrinsicDelaunayTriangulation( const std::shared_ptr<const CombinatorialMap>& base,
                                        const VertexPositionsFunc& vert_positions );

        virtual ~IntrinsicDelaunayTriangulation() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override { return mBaseMap->maxDartId(); }

        virtual uint dim() const override { return 2; }

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        /// @brief Returns the 3D position of the vertex represented by a given dart in this map.
        Eigen::Vector3d vertexPosition( const Vertex& v ) const;
        double edgeLength( const Edge& e ) const;
        double cotangentWeight( const Edge& e ) const;

        const CombinatorialMap& baseMap() const { return *mBaseMap; }

        private:
        const std::shared_ptr<const CombinatorialMap> mBaseMap;
        std::unordered_map<Dart::IndexType, Dart> mAlteredPhi1s;
        std::unordered_map<Dart::IndexType, Dart> mAlteredPhi_1s;

        /// Maps each dart id to its current vertex id in the intrinsic Delaunay map.
        std::unordered_map<Dart::IndexType, size_t> mDartToVertex;

        /// Physical 3D positions of vertices, indexed by vertex id from the base map's indexing.
        std::vector<Eigen::Vector3d> mVertexPositions;
        std::unordered_map<Dart::IndexType, size_t> mEdgeIds;
        std::vector<double> mEdgeLengths;

        Dart phi1Current( const Dart& d ) const;
        Dart phi_1Current( const Dart& d ) const;

        void setPhi1( const Dart& d, const Dart& target );
        void flipEdge( const Dart& d );
        bool isDelaunay( const Dart& d ) const;
        size_t edgeId( const Dart& d ) const;
        double edgeLength( const Dart& d ) const;
    };

    /// @brief Returns a vertex positions function for an IntrinsicDelaunayTriangulation that
    /// correctly resolves the 3D position for each vertex in the modified topology.
    VertexPositionsFunc intrinsicDelaunayVertexPositions( const IntrinsicDelaunayTriangulation& idtri );

} // namespace topology
