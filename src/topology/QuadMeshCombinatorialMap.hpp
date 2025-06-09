#pragma once
#include <SimplicialComplex.hpp>
#include <map>
#include <CombinatorialMap.hpp>

namespace topology
{
    /// @brief A combinatorial map describing an all-Quad mesh.
    /// The phi operations are defined by using a canonical Quad-local
    /// dart ordering, so that only the phi2 operations need O(n)
    /// data.
    /// See "Compact Combinatorial Maps: a Volume Mesh Data Structure,"
    /// by Xin Feng, et al., which describes this technique. Instead of
    /// storing a lookup table for phi1s as described in that paper, we
    /// set the local dart ordering such that phi1s may be calculated
    /// using simple arithmetic.
    ///
    /// Quad local dart structure, relative to local vertex ids:
    ///  2 ________________ 1
    ///   |                |
    ///   | o  ---------o  |
    ///   | |      1       |
    ///   | |           |  |
    ///   | | 2         |  |
    ///   | |         0 |  |
    ///   | |           |  |
    ///   |       3     |  |
    ///   | o---------  o  |
    ///   |________________|
    ///  3                  0
    class QuadMeshCombinatorialMap : public CombinatorialMap
    {
        public:
        QuadMeshCombinatorialMap( const std::vector<std::array<VertexId, 4>>& quads, const size_t num_vertices );

        virtual ~QuadMeshCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override;
        virtual uint dim() const override { return 2; }
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        Vertex vertexOfId( const VertexId& vid ) const
        {
            return Vertex( mVertexDarts.at( vid.id() ) );
        }

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        private:
        Dart phi1_1( const int i, const Dart& ) const;
        const std::vector<std::array<VertexId, 4>> mVerticesOfQuads;

        /// Stores one phi2 for each half edge.  The other phi2s can be
        /// calculated using these and phi1 operations.  The mth entry of
        /// this vector is the phi2 of the dart with id 4*m.
        std::vector<Dart> mPhi2s;

        std::vector<size_t> mEdgeIds;
        std::vector<Dart> mVertexDarts;

        size_t mNumEdges;

        Dart dartOfQuad( const uint quad_id ) const;
    };

} // namespace topology