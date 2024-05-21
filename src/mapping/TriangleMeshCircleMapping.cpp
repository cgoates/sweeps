#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <numbers>

namespace mapping
{
    TriangleMeshCircleMapping::TriangleMeshCircleMapping( const param::TriangleParametricAtlas& atlas, const VertexPositionsFunc& vertex_positions ) :
        mAtlas( atlas ),
        mPositions( vertex_positions )
    {
        const auto vert_ids = indexingOrError( mAtlas.cmap(), 0 );
        iterateCellsWhile( mAtlas.cmap(), 0, [&]( const topology::Vertex& v ) {
            if( boundaryAdjacent( mAtlas.cmap(), v ) )
            {
                const Eigen::Vector2d pos = mPositions( v );
                mBoundaryAngles.emplace( vert_ids( v ), atan2( pos( 1 ), pos( 0 ) ) );
            }
            return true;
        } );
    }

    std::optional<topology::Edge> maybeBoundaryEdge( const topology::CombinatorialMap& cmap, const topology::Face& f )
    {
        std::optional<topology::Edge> out;
        iterateDartsOfCell( cmap, f, [&]( const topology::Dart& d ){
            if( not phi( cmap, cmap.dim(), d ).has_value() )
            {
                out.emplace( d );
                return false;
            }
            return true;
        } );
        return out;
    }

    Eigen::VectorXd TriangleMeshCircleMapping::evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const
    {
        if( c.dim() != 2 ) throw std::runtime_error( "Bad cell dimension for TriangleMeshCircleMapping::evaluate" );
        const auto vertex_ii = [&]( const topology::Vertex& v ) {
            const param::ParentPoint v_pt = mAtlas.parentPoint( v );
            return std::distance( v_pt.mBaryCoordIsZero.begin(),
                                  std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
        };

        const std::optional<topology::Edge> boundary_edge = maybeBoundaryEdge( mAtlas.cmap(), c );
        const Vector6dMax expanded_coords = expandedCoordinates( pt );
        if( boundary_edge.has_value() )
        {
            const auto vert_ids = indexingOrError( mAtlas.cmap(), 0 );
            const topology::Vertex v0( boundary_edge.value().dart() );
            const topology::Vertex v1( phi( mAtlas.cmap(), 1, boundary_edge.value().dart() ).value() );
            const topology::Vertex v2( phi( mAtlas.cmap(), -1, boundary_edge.value().dart() ).value() );
            const double& bary0 = expanded_coords( vertex_ii( v0 ) );
            const double& bary1 = expanded_coords( vertex_ii( v1 ) );
            const double& bary2 = expanded_coords( vertex_ii( v2 ) );
            // 1. Get the theta value
            const double theta =
                ( mBoundaryAngles.at( vert_ids( v0 ) ) * bary0 + mBoundaryAngles.at( vert_ids( v1 ) ) * bary1 ) /
                ( bary0 + bary1 );

            // Check for bary0 + bary1 = 0 without toleranced compares.  Bad idea?
            if( std::isinf( theta ) or std::isnan( theta ) ) return mPositions( v2 );

            // 2. Evaluate the arc at that value
            const Eigen::Vector2d boundary_point( cos( theta ), sin( theta ) );

            // 3. Interpolate between the arc and the third point
            return ( bary0 + bary1 ) * boundary_point + bary2 * mPositions( v2 );
        }
        else
        {
            Vector3dMax out = Vector3dMax::Zero( 2 );
            iterateAdjacentCellsOfRestrictedCell( mAtlas.cmap(), c, 2, 0, [&]( const topology::Vertex& v ) {
                out += mPositions( v ) * expanded_coords( vertex_ii( v ) );
                return true;
            } );
            return out;
        }
    }

}