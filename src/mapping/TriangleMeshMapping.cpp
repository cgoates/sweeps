#include <TriangleMeshMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SimplexUtilities.hpp>
#include <CommonUtils.hpp>

namespace mapping
{
    TriangleMeshMapping::TriangleMeshMapping( const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
                                              const VertexPositionsFunc& vertex_positions,
                                              const size_t dim )
        : mAtlas( atlas ), mPositions( vertex_positions ), mDim( dim )
    {}

    Eigen::VectorXd TriangleMeshMapping::evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const
    {
        if( c.dim() != 2 ) throw std::runtime_error( "Bad cell dimension for TriangleMeshMapping::evaluate" );
        const auto vertex_ii = [&]( const topology::Vertex& v ) {
            const param::ParentPoint v_pt = mAtlas->parentPoint( v );
            return std::distance( v_pt.mBaryCoordIsZero.begin(),
                                  std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
        };

        const Vector6dMax expanded_coords = expandedCoordinates( pt );
        Vector3dMax out;
        bool first = true;
        iterateAdjacentCellsOfRestrictedCell( mAtlas->cmap(), c, 2, 0, [&]( const topology::Vertex& v ) {
            if( first )
            {
                out = mPositions( v ) * expanded_coords( vertex_ii( v ) );
                first = false;
            }
            else out += mPositions( v ) * expanded_coords( vertex_ii( v ) );
            return true;
        } );
        return out;
    }

    std::optional<param::ParentPoint> TriangleMeshMapping::maybeInverse( const topology::Face& f, const Eigen::Vector2d& pt ) const
    {
        const auto vertex_ii = [&]( const topology::Vertex& v ) {
            const param::ParentPoint v_pt = mAtlas->parentPoint( v );
            return std::distance( v_pt.mBaryCoordIsZero.begin(),
                                  std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
        };
        
        std::optional<param::ParentPoint> out = std::nullopt;

        iterateAdjacentCellsOfRestrictedCell( mAtlas->cmap(), f, 2, 0, [&]( const topology::Vertex& v ) {
            if( vertex_ii( v ) == 0 )
            {
                const Triangle<2> tri = triangleOfFace<2>( mAtlas->cmap(), mPositions, topology::Face( v.dart() ) );
                const param::ParentDomain pd = mAtlas->parentDomain( f );
                const std::optional<Eigen::Vector3d> maybe_bary = invertTriangleMap( tri, pt );
                out = maybe_bary.transform( [&]( const Eigen::Vector3d& bary_coords ) {
                    return compressCoordinates( pd, bary_coords, 1e-5 );
                } );

                return false;
            }
            return true;
        } );
        return out;
    }

    std::optional<std::pair<topology::Face, param::ParentPoint>> TriangleMeshMapping::maybeInverse( const Eigen::Vector2d& pt ) const
    {
        std::optional<std::pair<topology::Face, param::ParentPoint>> out;
        iterateCellsWhile( mAtlas->cmap(), 2, [&]( const topology::Face& f ) {
            out = maybeInverse( f, pt ).transform( [&f]( const param::ParentPoint& ppt ) {
                return std::pair<topology::Face, param::ParentPoint>{ f, ppt };
            } );
            return not out.has_value();
        } );
        return out;
    }

}