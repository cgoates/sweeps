#pragma once
#include <SimplicialComplex.hpp>
#include <map>
#include <CombinatorialMap.hpp>

namespace topology
{
    /// @brief A combinatorial map describing an all-tri mesh.
    /// The phi operations are defined by using a canonical tri-local
    /// dart ordering, so that only the phi2 operations need O(n)
    /// data.
    /// See "Compact Combinatorial Maps: a Volume Mesh Data Structure,"
    /// by Xin Feng, et al., which describes this technique. Instead of
    /// storing a lookup table for phi1s as described in that paper, we
    /// set the local dart ordering such that phi1s may be calculated
    /// using simple arithmetic.
    ///
    /// Tri local dart structure, relative to local vertex ids:
    ///                    1
    ///                   /\.
    ///                  /  \.
    ///                 / o  \.
    ///                / /    \.
    ///               / / 1    \.
    ///              / /      \ \.
    ///             /    2   0 \ \.
    ///            /  o----     o \.
    ///         2 /________________\ 0
    ///
    class TriMeshCombinatorialMap : public CombinatorialMap
    {
        public:
        TriMeshCombinatorialMap( const SimplicialComplex& complex );

        virtual ~TriMeshCombinatorialMap() = default;

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

        const SimplicialComplex& simplicialComplex() const { return mSimplicialComplex; }

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        private:
        const SimplicialComplex& mSimplicialComplex;

        /// Stores one phi3 for each half face.  The other phi3s can be
        /// calculated using these and phi1 operations.  The mth entry of
        /// this vector is the phi3 of the dart with id 3*m.
        std::vector<Dart> mPhi2s;

        std::vector<size_t> mEdgeIds;
        std::vector<Dart> mVertexDarts;

        size_t mNumEdges;

        Dart dartOfTri( const uint tri_id ) const;
    };

} // namespace topology