#include <IntrinsicDelaunayTriangulation.hpp>
#include <CombinatorialMapMethods.hpp>
#include <GlobalCellMarker.hpp>
#include <queue>
#include <cmath>
#include <numbers>

using namespace topology;

static double clampedCosineFromLengths( const double adj1, const double adj2, const double opp )
{
    const double denom = 2.0 * adj1 * adj2;
    if( denom <= 0.0 ) throw std::runtime_error( "Degenerate intrinsic triangle" );
    return std::clamp( ( adj1 * adj1 + adj2 * adj2 - opp * opp ) / denom, -1.0, 1.0 );
}

static double angleFromLengths( const double adj1, const double adj2, const double opp )
{
    return std::acos( clampedCosineFromLengths( adj1, adj2, opp ) );
}

static double cotangentFromLengths( const double adj1, const double adj2, const double opp )
{
    const double cos_theta = clampedCosineFromLengths( adj1, adj2, opp );
    const double sin_sq = std::max( 1.0 - cos_theta * cos_theta, 1e-30 );
    return cos_theta / std::sqrt( sin_sq );
}

static double edgeLengthFromPositions( const Eigen::Vector3d& p1, const Eigen::Vector3d& p2 )
{
    return ( p2 - p1 ).norm();
}

Dart IntrinsicDelaunayTriangulation::phi1Current( const Dart& d ) const
{
    const auto it = mAlteredPhi1s.find( d.id() );
    if( it != mAlteredPhi1s.end() ) return it->second;
    return topology::phi( *mBaseMap, 1, d ).value();
}

Dart IntrinsicDelaunayTriangulation::phi_1Current( const Dart& d ) const
{
    const auto it = mAlteredPhi_1s.find( d.id() );
    if( it != mAlteredPhi_1s.end() ) return it->second;
    return topology::phi( *mBaseMap, -1, d ).value();
}

void IntrinsicDelaunayTriangulation::setPhi1( const Dart& d, const Dart& target )
{
    mAlteredPhi1s[d.id()] = target;
    mAlteredPhi_1s[target.id()] = d;
}

void IntrinsicDelaunayTriangulation::flipEdge( const Dart& a )
{
    // Collect the 6 darts of the two triangles sharing edge a before modifying anything.
    const Dart b       = phi1Current( a );               // next in F1
    const Dart c       = phi_1Current( a );              // prev in F1 (= opposite vertex k)
    const Dart a_prime = topology::phi( *mBaseMap, 2, a ).value(); // phi2, always from base
    const Dart b_prime = phi1Current( a_prime );         // next in F2
    const Dart c_prime = phi_1Current( a_prime );        // prev in F2 (= opposite vertex l)

    // Update the intrinsic length on the flipped diagonal by unfolding the two incident
    // triangles into a common plane around vertex i.
    const double l_ij = edgeLength( a );
    const double l_jk = edgeLength( b );
    const double l_ki = edgeLength( c );
    const double l_il = edgeLength( b_prime );
    const double l_lj = edgeLength( c_prime );
    const double angle_kil = angleFromLengths( l_ij, l_ki, l_jk ) + angleFromLengths( l_ij, l_il, l_lj );
    const double new_len_sq =
        std::max( 0.0, l_ki * l_ki + l_il * l_il - 2.0 * l_ki * l_il * std::cos( angle_kil ) );
    const double new_len = std::sqrt( new_len_sq );
    mEdgeLengths.at( edgeId( a ) ) = new_len;

    // Update vertex assignments: only a and a_prime change.
    // After flip: Vertex(a) = Vertex(c), Vertex(a') = Vertex(c').
    mDartToVertex[a.id()]       = mDartToVertex.at( c.id() );
    mDartToVertex[a_prime.id()] = mDartToVertex.at( c_prime.id() );

    // Rewire phi1:
    //   New face cycle 1: a → c' → b → a
    //   New face cycle 2: a' → c → b' → a'
    setPhi1( a,       c_prime );
    setPhi1( c_prime, b       );
    setPhi1( b,       a       );
    setPhi1( a_prime, c       );
    setPhi1( c,       b_prime );
    setPhi1( b_prime, a_prime );
}

bool IntrinsicDelaunayTriangulation::isDelaunay( const Dart& d ) const
{
    const auto maybe_adj = topology::phi( *mBaseMap, 2, d );
    if( !maybe_adj.has_value() ) return true; // boundary edge is always Delaunay

    const Dart a_prime = maybe_adj.value();
    const Dart b = phi1Current( d );
    const Dart c = phi_1Current( d );
    const Dart b_prime = phi1Current( a_prime );
    const Dart c_prime = phi_1Current( a_prime );

    const double l_ij = edgeLength( d );
    const double l_jk = edgeLength( b );
    const double l_ki = edgeLength( c );
    const double l_il = edgeLength( b_prime );
    const double l_lj = edgeLength( c_prime );

    const double alpha = angleFromLengths( l_ki, l_jk, l_ij );
    const double beta = angleFromLengths( l_il, l_lj, l_ij );
    constexpr double eps = 1e-10;
    return alpha + beta <= std::numbers::pi + eps;
}

IntrinsicDelaunayTriangulation::IntrinsicDelaunayTriangulation(
    const std::shared_ptr<const CombinatorialMap>& base,
    const VertexPositionsFunc& vert_positions )
    : mBaseMap( base )
{
    // Initialise mDartToVertex and mVertexPositions from the base map.
    const auto maybe_base_indexing = mBaseMap->indexing( 0 );
    if( !maybe_base_indexing.has_value() )
        throw std::runtime_error( "IntrinsicDelaunayTriangulation: base map must provide vertex indexing" );
    const auto base_indexing = maybe_base_indexing.value();

    const size_t n_verts = topology::cellCount( *mBaseMap, 0 );
    mVertexPositions.resize( n_verts );
    mEdgeLengths.resize( topology::cellCount( *mBaseMap, 1 ) );

    topology::iterateCellsWhile( *mBaseMap, 0, [&]( const Vertex& v ) {
        const size_t vid = base_indexing( v );
        mVertexPositions[vid] = vert_positions( v );
        iterateDartsOfCell( *mBaseMap, v, [&]( const Dart& d ) {
            mDartToVertex[d.id()] = vid;
            return true;
        } );
        return true;
    } );
    topology::iterateCellsWhile( *mBaseMap, 1, [&]( const Edge& e ) {
        mEdgeIds.emplace( lowestDartId( *mBaseMap, e ), mEdgeIds.size() );
        const Dart d = e.dart();
        const double len = edgeLengthFromPositions(
            mVertexPositions.at( mDartToVertex.at( d.id() ) ),
            mVertexPositions.at( mDartToVertex.at( phi1Current( d ).id() ) ) );
        mEdgeLengths.at( edgeId( d ) ) = len;
        return true;
    } );

    // Iterative edge-flip loop using a queue.
    std::queue<Dart::IndexType> q;
    std::vector<bool> queued_edges( mEdgeLengths.size(), false );
    const auto queue_edge = [&]( const Dart& d ) {
        if( !topology::phi( *mBaseMap, 2, d ).has_value() ) return;
        const size_t eid = edgeId( d );
        if( !queued_edges.at( eid ) )
        {
            queued_edges[eid] = true;
            q.push( d.id() );
        }
    };
    topology::iterateCellsWhile( *mBaseMap, 1, [&]( const Edge& e ) {
        queue_edge( e.dart() );
        return true;
    } );

    while( !q.empty() )
    {
        const Dart d( q.front() );
        q.pop();
        queued_edges.at( edgeId( d ) ) = false;

        // Re-check: edge might have become Delaunay since it was queued.
        if( isDelaunay( d ) ) continue;

        // Collect adjacent darts BEFORE the flip.
        const Dart b       = phi1Current( d );
        const Dart c       = phi_1Current( d );
        const Dart a_prime = topology::phi( *mBaseMap, 2, d ).value();
        const Dart b_prime = phi1Current( a_prime );
        const Dart c_prime = phi_1Current( a_prime );

        flipEdge( d );

        // Add the 4 surrounding interior edges to the queue.
        for( const Dart& adj : { b, c, b_prime, c_prime } )
        {
            queue_edge( adj );
        }
    }
}

std::optional<Dart> IntrinsicDelaunayTriangulation::phi( const int i, const Dart& d ) const
{
    if( i == 1 )
    {
        const auto it = mAlteredPhi1s.find( d.id() );
        if( it != mAlteredPhi1s.end() ) return it->second;
        return topology::phi( *mBaseMap, 1, d );
    }
    else if( i == -1 )
    {
        const auto it = mAlteredPhi_1s.find( d.id() );
        if( it != mAlteredPhi_1s.end() ) return it->second;
        return topology::phi( *mBaseMap, -1, d );
    }
    else if( i == 2 )
    {
        return topology::phi( *mBaseMap, 2, d );
    }
    throw std::runtime_error( "IntrinsicDelaunayTriangulation: bad phi index" );
}

bool IntrinsicDelaunayTriangulation::iterateDartsWhile(
    const std::function<bool( const Dart& )>& callback ) const
{
    return topology::iterateDartsWhile( *mBaseMap, callback );
}

bool IntrinsicDelaunayTriangulation::iterateCellsWhile(
    const uint cell_dim,
    const std::function<bool( const Cell& )>& callback ) const
{
    const auto indexing = this->indexing( cell_dim );
    if( indexing.has_value() )
    {
        IndexingCellMarker m( *indexing, cell_dim );
        return iterateDartsWhile( [&]( const Dart& d ) {
            const Cell c( d, cell_dim );
            if( !m.isMarked( c ) )
            {
                m.mark( c );
                if( !callback( c ) ) return false;
            }
            return true;
        } );
    }

    GlobalCellMarker m( *this, cell_dim );
    return iterateDartsWhile( [&]( const Dart& d ) {
        const Cell c( d, cell_dim );
        if( !m.isMarked( c ) )
        {
            m.mark( *this, c );
            if( !callback( c ) ) return false;
        }
        return true;
    } );
}

std::optional<IndexingFunc> IntrinsicDelaunayTriangulation::indexing( const uint cell_dim ) const
{
    if( cell_dim == 0 )
    {
        return [this]( const Vertex& v ) -> size_t {
            const auto it = mDartToVertex.find( v.dart().id() );
            if( it != mDartToVertex.end() ) return it->second;
            throw std::runtime_error( "IntrinsicDelaunayTriangulation: dart not in vertex table" );
        };
    }
    else if( cell_dim == 1 || cell_dim == 2 )
    {
        return mBaseMap->indexing( cell_dim );
    }
    return std::nullopt;
}

std::optional<size_t> IntrinsicDelaunayTriangulation::cellCount( const uint cell_dim ) const
{
    return mBaseMap->cellCount( cell_dim );
}

Eigen::Vector3d IntrinsicDelaunayTriangulation::vertexPosition( const Vertex& v ) const
{
    const auto it = mDartToVertex.find( v.dart().id() );
    if( it == mDartToVertex.end() )
        throw std::runtime_error( "IntrinsicDelaunayTriangulation::vertexPosition: dart not found" );
    return mVertexPositions.at( it->second );
}

double IntrinsicDelaunayTriangulation::edgeLength( const Dart& d ) const
{
    return mEdgeLengths.at( edgeId( d ) );
}

double IntrinsicDelaunayTriangulation::edgeLength( const Edge& e ) const
{
    return edgeLength( e.dart() );
}

size_t IntrinsicDelaunayTriangulation::edgeId( const Dart& d ) const
{
    return mEdgeIds.at( lowestDartId( *mBaseMap, Edge( d ) ) );
}

double IntrinsicDelaunayTriangulation::cotangentWeight( const Edge& e ) const
{
    const Dart d = e.dart();
    const double l_uv = edgeLength( d );
    const double cot1 = cotangentFromLengths( edgeLength( phi_1Current( d ) ), edgeLength( phi1Current( d ) ), l_uv );

    const auto maybe_adj = topology::phi( *mBaseMap, 2, d );
    if( !maybe_adj.has_value() ) return cot1;

    const Dart a_prime = maybe_adj.value();
    const double cot2 = cotangentFromLengths( edgeLength( phi1Current( a_prime ) ), edgeLength( phi_1Current( a_prime ) ), l_uv );
    return cot1 + cot2;
}

namespace topology
{
    VertexPositionsFunc intrinsicDelaunayVertexPositions( const IntrinsicDelaunayTriangulation& idtri )
    {
        const auto* idtri_ptr = &idtri;
        return [idtri_ptr]( const Vertex& v ) -> Eigen::Vector3d {
            return idtri_ptr->vertexPosition( v );
        };
    }
}
