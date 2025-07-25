#pragma once
#include <SimplicialComplex.hpp>
#include <map>
#include <CombinatorialMap.hpp>

namespace topology
{
    /// @brief A combinatorial map describing an all-tet mesh.
    /// The phi operations are defined by using a canonical tet-local
    /// dart ordering, so that only the phi3 operations need O(n)
    /// data.
    /// See "Compact Combinatorial Maps: a Volume Mesh Data Structure,"
    /// by Xin Feng, et al., which describes this technique. Instead of
    /// storing a lookup table for phi1s as described in that paper, we
    /// set the local dart ordering such that phi1s may be calculated
    /// using simple arithmetic.
    ///
    /// Tet local dart structure, relative to local vertex ids:
    ///                    3
    ///                   /\.
    ///                  /  \.
    ///                 / o  \.
    ///                / /    \.
    ///               / / 5    \.
    ///              / /      \ \.
    ///             /    3   4 \ \.
    ///            /  o----     o \.
    ///         1 /________________\ 0
    ///          /\ o     ------ o /\.
    ///         /  \ \ 1     0    /  \.
    ///        / o  \ \          / o  \.
    ///       / /    \        / / /    \.
    ///      / / 7  \ \    2 / / / 9  \ \.
    ///     / /      \ \    / / /      \ \.
    ///    /    8   6 \ \  o /   10  11 \ \.
    ///   /  o----     o \  /  o----     o \.
    ///  /________________\/________________\.
    /// 3                 2                  3
    ///
    /// The diagram above is to be read as if looking at the inside of
    /// the tet, i.e., all of the vertex 3 points should be folded out
    /// of the page to form the tetrahedron.
    class TetMeshCombinatorialMap : public CombinatorialMap
    {
        public:
        TetMeshCombinatorialMap( const SimplicialComplex& complex );

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;
        virtual Dart::IndexType maxDartId() const override;
        virtual uint dim() const override { return 3; }
        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;
        virtual bool iterateCellsWhile( const uint cell_dim, const std::function<bool( const Cell& )>& callback ) const override;
        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        Vertex vertexOfId( const VertexId& vid ) const
        {
            return Vertex( mVertexDarts.at( vid.id() ) );
        }

        Volume tetOfId( const size_t tet_id ) const;

        const SimplicialComplex& simplicialComplex() const { return mSimplicialComplex; }

        virtual std::optional<size_t> cellCount( const uint cell_dim ) const override;

        private:
        Dart phi1_1( const int, const Dart& ) const;
        const SimplicialComplex& mSimplicialComplex;

        /// This array holds the tet-local phi2 operations. The ith entry is
        /// the tet-local dart id of the phi2 of the dart with tet-local id i.
        /// e.g. phi2( 0 ) = 3, as seen in the diagram above.
        static constexpr std::array<Dart::IndexType, 12> mPhi2s = { 3, 6, 9, 0, 11, 7, 1, 5, 10, 2, 8, 4 };

        static constexpr std::array<VertexId::Type, 12> mLocalVertices = { 0, 1, 2, 1, 0, 3, 2, 1, 3, 0, 2, 3 };
        static constexpr std::array<Dart::IndexType, 4> mLocalFaceDarts = { 0, 3, 6, 9 };
        static constexpr std::array<Dart::IndexType, 4> mLocalVertexDarts = { 0, 1, 2, 5 };

        /// Stores one phi3 for each half face.  The other phi3s can be
        /// calculated using these and phi1 operations.  The mth entry of
        /// this vector is the phi3 of the dart with id 3*m.
        std::vector<Dart> mPhi3s;

        std::vector<size_t> mFaceIds;
        std::vector<size_t> mEdgeIds;
        std::vector<Dart> mVertexDarts;

        size_t mNumEdges;
        size_t mNumFaces;

        Dart dartOfTet( const uint tet_id ) const;
    };

} // namespace topology