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

        const Vector6dMax expanded_coords = expandedCoordinates( pt );
        Vector3dMax out;
        bool first = true;
        iterateAdjacentCellsOfRestrictedCell( mAtlas->cmap(), c, 2, 0, [&]( const topology::Vertex& v ) {
            if( first )
            {
                out = mPositions( v ) * expanded_coords( vertexIndex( v ) );
                first = false;
            }
            else out += mPositions( v ) * expanded_coords( vertexIndex( v ) );
            return true;
        } );
        return out;
    }

    std::optional<param::ParentPoint> TriangleMeshMapping::maybeInverse( const topology::Face& f, const Eigen::Vector2d& pt ) const
    {
        std::optional<param::ParentPoint> out = std::nullopt;

        iterateAdjacentCellsOfRestrictedCell( mAtlas->cmap(), f, 2, 0, [&]( const topology::Vertex& v ) {
            if( vertexIndex( v ) == 0 )
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

    std::optional<std::pair<topology::Cell, param::ParentPoint>> TriangleMeshMapping::maybeInverse( const Eigen::Vector2d& pt ) const
    {
        std::optional<std::pair<topology::Cell, param::ParentPoint>> out;
        iterateCellsWhile( mAtlas->cmap(), 2, [&]( const topology::Face& f ) {
            out = maybeInverse( f, pt ).transform( [&f]( const param::ParentPoint& ppt ) {
                return std::pair<topology::Cell, param::ParentPoint>{ f, ppt };
            } );
            return not out.has_value();
        } );
        return out;
    }

    std::pair<topology::Cell, param::ParentPoint> TriangleMeshMapping::closestPoint( const Eigen::Vector3d& pt ) const
    {
        std::optional<std::pair<topology::Cell, param::ParentPoint>> out;
        double min_dist = std::numeric_limits<double>::max();
        iterateCellsWhile( mAtlas->cmap(), 2, [&]( const topology::Face& f ) {
            bool on_triangle = false;
            iterateAdjacentCellsOfRestrictedCell( mAtlas->cmap(), f, 2, 0, [&]( const topology::Vertex& v ) {
                if( vertexIndex( v ) == 0 )
                {
                    const Triangle<3> tri = triangleOfFace<3>( mAtlas->cmap(), mPositions, topology::Face( v.dart() ) );
                    const std::pair<Eigen::Vector3d, std::optional<double>> maybe_closest = invertTriangleMapOrClosestPoint( tri, pt );
                    const param::ParentDomain pd = mAtlas->parentDomain( f );
                    if( not maybe_closest.second.has_value() )
                    {
                        out.emplace( std::pair<topology::Cell, param::ParentPoint>{ f, compressCoordinates( pd, maybe_closest.first, 1e-5 ) } );
                        on_triangle = true;
                    }
                    else
                    {
                        const double dist = maybe_closest.second.value();
                        if( dist < min_dist )
                        {
                            min_dist = dist;
                            out.emplace( std::pair<topology::Cell, param::ParentPoint>{ f, compressCoordinates( pd, maybe_closest.first, 1e-5 ) } );
                        }
                    }
                    return false;
                }
                return true;
            } );
            
            return not on_triangle;
        } );

        if( not out.has_value() ) throw std::runtime_error( "No closest point found" );

        return out.value();
    }

    size_t TriangleMeshMapping::vertexIndex( const topology::Vertex& v ) const
    {
        const param::ParentPoint v_pt = mAtlas->parentPoint( v );
        return std::distance( v_pt.mBaryCoordIsZero.begin(),
                              std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
    }
}