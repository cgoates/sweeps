#pragma once
#include <SimplicialComplex.hpp>
#include <map>

class Dart
{
    public:
    using IndexType = uint64_t;

    Dart( const IndexType& id ) : mIndex( id ) {}

    IndexType index() const { return mIndex; }

    bool operator<( const Dart& o ) const
    {
        return index() < o.index();
    }

    private:

    // TODO: Perhaps add a reference to the combinatorial map here
    // so we can check that darts are being passed to the correct map?
    // Maybe we don't need that. Hold off for now. It is nice to have
    // mutablity.

    IndexType mIndex;
};

class Cell
{
    public:
    Cell( const Dart& d, const uint dim ) : mDart( d ), mDim( dim ) {}

    const Dart& dart() const { return mDart; }
    uint dim() const { return mDim; }

    private:
    Dart mDart;
    uint mDim;
};

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
class TetMeshCombinatorialMap
{
    public:
    TetMeshCombinatorialMap( const SimplicialComplex& complex );

    std::optional<Dart> phi( const int i, const Dart& d ) const;

    size_t elementId( const Cell& c ) const;
    VertexId vertexId( const Cell& c ) const;

    private:
    const SimplicialComplex& mSimplicialComplex;

    /// This array holds the tet-local phi2 operations. The ith entry is
    /// the tet-local dart id of the phi2 of the dart with tet-local id i.
    /// e.g. phi2( 0 ) = 3, as seen in the diagram above.
    static constexpr std::array<Dart::IndexType, 12> mPhi2s = { 3, 6, 9, 0, 11, 7, 1, 5, 10, 2, 8, 4 };

    /// Stores one phi3 for each half face.  The other phi3s can be
    /// calculated using these and phi1 operations.
    std::map<Dart, Dart> mPhi3s;

    static constexpr std::array<VertexId::Type, 12> mLocalVertices = { 0, 1, 2, 1, 0, 3, 2, 1, 3, 0, 2, 3 };
};