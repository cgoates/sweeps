#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <numbers>
#include <SimplexUtilities.hpp>

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
            const double theta0 = mBoundaryAngles.at( vert_ids( v0 ) );
            const double theta1 = [&]() {
                double temp = mBoundaryAngles.at( vert_ids( v1 ) );
                while( temp < theta0 ) temp += 2 * std::numbers::pi;
                return temp;
            }();
            // 1. Get the theta value
            const double theta = ( theta0 * bary0 + theta1 * bary1 ) / ( bary0 + bary1 );

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

    // Intersection of the line formed by pt0 and pt1 with the unit circle,
    // in the direction of pt1 - pt0. Assumes that both points are in the circle.
    Eigen::Vector2d lineCircleIntersection( const Eigen::Vector2d& pt0, const Eigen::Vector2d& pt1 )
    {
        const double A = pt1( 1 ) - pt0( 1 );
        const double B = pt0( 0 ) - pt1( 0 );
        const double C = -A * pt0( 0 ) - B * pt0( 1 );

        const double Asquared = A*A;
        const double Bsquared = B*B;
        const double Csquared = C*C;

        const double denominator = 1.0 / ( Asquared + Bsquared );
        const double sqrt_term = sqrt( Asquared + Bsquared - Csquared );

        return Eigen::Vector2d( -( sqrt_term * B + A * C ) * denominator, ( sqrt_term * A - B * C ) * denominator );
    }

    double normalizeAngle( double angle )
    {
        while( angle < -1 * std::numbers::pi ) angle += 2 * std::numbers::pi;
        while( angle > std::numbers::pi ) angle -= 2 * std::numbers::pi;
        return angle;
    }

    // See https://stackoverflow.com/a/23550032
    bool isBetweenAngles( double a, double b, const double test_angle )
    {
        a -= test_angle;
        b -= test_angle;
        a = normalizeAngle( a );
        b = normalizeAngle( b );
        if( a * b > 0 ) return false;
        return std::abs( a - b ) < std::numbers::pi;
    }

    std::optional<param::ParentPoint> TriangleMeshCircleMapping::maybeInverse( const topology::Face& f, const Eigen::Vector2d& pt ) const
    {
        const auto vertex_ii = [&]( const topology::Vertex& v ) {
            const param::ParentPoint v_pt = mAtlas.parentPoint( v );
            return std::distance( v_pt.mBaryCoordIsZero.begin(),
                                  std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
        };
        const std::optional<topology::Edge> boundary_edge = maybeBoundaryEdge( mAtlas.cmap(), f );
        if( boundary_edge.has_value() )
        {
            const topology::Vertex v0( boundary_edge.value().dart() );
            const topology::Vertex v1( phi( mAtlas.cmap(), 1, v0.dart() ).value() );
            const topology::Vertex v2( phi( mAtlas.cmap(), 1, v1.dart() ).value() );
            const Eigen::Vector2d non_bdry_point = mPositions( v2 );

            const Eigen::Vector2d one_bdry_vert_point = mPositions( v0 );
            if( ( pt - non_bdry_point ).norm() < 1e-10 * ( one_bdry_vert_point - non_bdry_point ).norm() ) return mAtlas.parentPoint( v2 );

            const auto vert_ids = indexingOrError( mAtlas.cmap(), 0 );
            const double theta0 = mBoundaryAngles.at( vert_ids( v0 ) );
            const double theta1 = mBoundaryAngles.at( vert_ids( v1 ) );

            const Eigen::Vector2d bdry_point = lineCircleIntersection( non_bdry_point, pt );
            const double bdry_point_theta = atan2( bdry_point( 1 ), bdry_point( 0 ) );

            if( not isBetweenAngles( theta0, theta1, bdry_point_theta ) ) return std::nullopt;

            const auto [ lambda_non_bdry, lambda_bdry ] = [&](){
                if( std::abs( bdry_point( 0 ) - non_bdry_point( 0 ) ) > std::abs( bdry_point( 1 ) - non_bdry_point( 1 ) ) )
                    return inverseLinear( non_bdry_point(0), bdry_point(0), pt(0) );
                else
                    return inverseLinear( non_bdry_point(1), bdry_point(1), pt(1) );
            }();

            const auto [lambda_v0, lambda_v1] = inverseLinear(
                normalizeAngle( theta0 - bdry_point_theta ), normalizeAngle( theta1 - bdry_point_theta ), 0.0 );

            Eigen::Vector3d bary = Eigen::Vector3d::Zero();
            bary( vertex_ii( v0 ) ) = lambda_v0 * lambda_bdry;
            bary( vertex_ii( v1 ) ) = lambda_v1 * lambda_bdry;
            bary( vertex_ii( v2 ) ) = lambda_non_bdry;

            const param::ParentDomain pd = mAtlas.parentDomain( f );
            return compressCoordinates( pd, bary, 1e-5 );
        }
        else
        {
            std::optional<param::ParentPoint> out = std::nullopt;

            iterateAdjacentCellsOfRestrictedCell( mAtlas.cmap(), f, 2, 0, [&]( const topology::Vertex& v ) {
                if( vertex_ii( v ) == 0 )
                {
                    const Triangle<2> tri = triangleOfFace<2>( mAtlas.cmap(), mPositions, topology::Face( v.dart() ) );
                    const param::ParentDomain pd = mAtlas.parentDomain( f );
                    const std::optional<Eigen::Vector3d> maybe_bary = invertTriangleMap( tri, pt );
                    out = maybe_bary.and_then( [&]( const Eigen::Vector3d& bary_coords ) -> std::optional<param::ParentPoint> {
                        return compressCoordinates( pd, bary_coords, 1e-5 );
                    } );

                    return false;
                }
                return true;
            } );
            return out;
        }
    }

    std::optional<std::pair<topology::Face, param::ParentPoint>> TriangleMeshCircleMapping::maybeInverse( const Eigen::Vector2d& pt ) const
    {
        if( pt.norm() > 1.0 + 1e-15 ) return std::nullopt;
        std::optional<std::pair<topology::Face, param::ParentPoint>> out;
        iterateCellsWhile( mAtlas.cmap(), 2, [&]( const topology::Face& f ) {
            out = maybeInverse( f, pt ).and_then(
                [&f]( const param::ParentPoint& ppt ) -> std::optional<std::pair<topology::Face, param::ParentPoint>> {
                    return std::pair<topology::Face, param::ParentPoint>{ f, ppt };
                } );
            return not out.has_value();
        } );
        return out;
    }

}