#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SimplexUtilities.hpp>
#include <CommonUtils.hpp>
#include <iostream>

namespace mapping
{
    SmallVector<topology::Edge, 3> maybeBoundaryEdges( const topology::CombinatorialMap& cmap, const topology::Face& f )
    {
        SmallVector<topology::Edge, 3> out;
        topology::Dart curr_d = f.dart();

        // Always start with a non-boundary dart so the boundary darts will be grouped together.
        while( onBoundary( cmap, curr_d ) ) curr_d = phi( cmap, -1, curr_d ).value();

        const topology::Dart ref_dart = curr_d;

        do
        {
            if( onBoundary( cmap, curr_d ) ) out.push_back( curr_d );

            curr_d = phi( cmap, 1, curr_d ).value();
        } while( curr_d != ref_dart );

        return out;
    }

    TriangleMeshCircleMapping::TriangleMeshCircleMapping(
        const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
        const VertexPositionsFunc& vertex_positions )
        : mAtlas( atlas ), mTriMapping( atlas, vertex_positions, 2, false )
    {
        const auto vert_ids = indexingOrError( mAtlas->cmap(), 0 );
        iterateCellsWhile( mAtlas->cmap(), 0, [&]( const topology::Vertex& v ) {
            if( boundaryAdjacent( mAtlas->cmap(), v ) )
            {
                const Eigen::Vector2d pos = mTriMapping.vertPositions()( v );
                mBoundaryAngles.emplace( vert_ids( v ), atan2( pos( 1 ), pos( 0 ) ) );
            }
            return true;
        } );

        iterateCellsWhile( mAtlas->cmap(), 2, [&]( const topology::Face& f ) {
            const Triangle<2> tri = triangleOfFace<2>( mAtlas->cmap(), vertex_positions, f );
            Eigen::Vector2d mins = tri.v1.array().min( tri.v2.array() ).min( tri.v3.array() );
            Eigen::Vector2d maxs = tri.v1.array().max( tri.v2.array() ).max( tri.v3.array() );

            const SmallVector<topology::Edge, 3> boundary_edges = maybeBoundaryEdges( mAtlas->cmap(), f );

            if( not boundary_edges.empty() )
            {
                const double max_angle = mBoundaryAngles.at( vert_ids( topology::Vertex( boundary_edges.front().dart() ) ) );
                const double min_angle = mBoundaryAngles.at(
                    vert_ids( topology::Vertex( phi( mAtlas->cmap(), 1, boundary_edges.back().dart() ).value() ) ) );

                // NOTE: We are assuming that the topological normal is in the same direction as the circle normal by RHR.
                // Nowhere do I check this, but that is how it is used so far. (4/2025)

                if( min_angle > max_angle )
                {
                    mins = mins.array().min( Eigen::Vector2d( -1, 0 ).array() );
                    maxs = maxs.array().max( Eigen::Vector2d( -1, 0 ).array() );
                }
                if( min_angle < 0 and max_angle > 0 )
                {
                    mins = mins.array().min( Eigen::Vector2d( 1, 0 ).array() );
                    maxs = maxs.array().max( Eigen::Vector2d( 1, 0 ).array() );
                }
                const double pi_half = std::numbers::pi / 2;
                if( min_angle < pi_half and max_angle > pi_half )
                {
                    mins = mins.array().min( Eigen::Vector2d( 0, 1 ).array() );
                    maxs = maxs.array().max( Eigen::Vector2d( 0, 1 ).array() );
                }
                if( min_angle < -pi_half and max_angle > -pi_half )
                {
                    mins = mins.array().min( Eigen::Vector2d( 0, -1 ).array() );
                    maxs = maxs.array().max( Eigen::Vector2d( 0, -1 ).array() );
                }
            }

            mBoundingBoxes.emplace( f, AABB( mins - Eigen::Vector2d::Constant( 1e-15 ), maxs + Eigen::Vector2d::Constant( 1e-15 ) ) );
            return true;
        } );
    }

    size_t vertex_ii( const param::ParametricAtlas& atlas, const topology::Vertex& v )
    {
        const param::ParentPoint v_pt = atlas.parentPoint( v );
        return std::distance( v_pt.mBaryCoordIsZero.begin(),
                              std::find( v_pt.mBaryCoordIsZero.begin(), v_pt.mBaryCoordIsZero.end(), false ) );
    }

    Eigen::Vector2d triangleWithArcEdgeEval( const std::array<double, 3>& bary, const std::array<double, 2>& thetas, const Eigen::Vector2d& v2_position )
    {
        const double& theta0 = thetas[0];
        const double theta1 = [&]() {
            double temp = thetas[1];
            while( temp < theta0 ) temp += 2 * std::numbers::pi;
            return temp;
        }();
        // 1. Get the theta value
        const double theta = ( theta0 * bary[0] + theta1 * bary[1] ) / ( bary[0] + bary[1] );

        if( std::isinf( theta ) or std::isnan( theta ) ) throw std::runtime_error( "Bad parent point provided to TriangleMeshCircleMapping::evaluate" );

        // 2. Evaluate the arc at that value
        const Eigen::Vector2d boundary_point( cos( theta ), sin( theta ) );

        // 3. Interpolate between the arc and the third point
        return ( bary[0] + bary[1] ) * boundary_point + bary[2] * v2_position;
    }

    Eigen::VectorXd TriangleMeshCircleMapping::evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const
    {
        if( c.dim() != 2 ) throw std::runtime_error( "Bad cell dimension for TriangleMeshCircleMapping::evaluate" );

        const SmallVector<topology::Edge, 3> boundary_edges = maybeBoundaryEdges( mAtlas->cmap(), c );
        switch( boundary_edges.size() )
        {
            case 0: return mTriMapping.evaluate( c, pt );
            case 1:
            {
                const topology::Edge& boundary_edge = boundary_edges.at( 0 );
                const Vector6dMax expanded_coords = expandedCoordinates( pt );
                const auto vert_ids = indexingOrError( mAtlas->cmap(), 0 );
                const std::array<topology::Vertex, 3> vertices{
                    topology::Vertex( boundary_edge.dart() ),
                    topology::Vertex( phi( mAtlas->cmap(), 1, boundary_edge.dart() ).value() ),
                    topology::Vertex( phi( mAtlas->cmap(), -1, boundary_edge.dart() ).value() )
                };
                const std::array<size_t, 3> corner_ids{
                    vertex_ii( *mAtlas, vertices[0] ),
                    vertex_ii( *mAtlas, vertices[1] ),
                    vertex_ii( *mAtlas, vertices[2] )
                };

                const std::array<double, 3> bary{
                    expanded_coords( vertex_ii( *mAtlas, vertices[0] ) ),
                    expanded_coords( vertex_ii( *mAtlas, vertices[1] ) ),
                    expanded_coords( vertex_ii( *mAtlas, vertices[2] ) )
                };
                const double theta0 = mBoundaryAngles.at( vert_ids( vertices[0] ) );
                const double theta1 = mBoundaryAngles.at( vert_ids( vertices[1] ) );

                if( pt.mBaryCoordIsZero.at( corner_ids[0] ) and pt.mBaryCoordIsZero.at( corner_ids[1] ) )
                    return mTriMapping.vertPositions()( vertices[2] );

                return triangleWithArcEdgeEval( bary, {theta0, theta1},  mTriMapping.vertPositions()( vertices[2] ) );
            }
            case 2:
            {
                const Vector6dMax expanded_coords = expandedCoordinates( pt );
                const auto vert_ids = indexingOrError( mAtlas->cmap(), 0 );
                const std::array<topology::Vertex, 3> vertices{
                    topology::Vertex( boundary_edges.at( 0 ).dart() ),
                    topology::Vertex( boundary_edges.at( 1 ).dart() ),
                    topology::Vertex( phi( mAtlas->cmap(), 1, boundary_edges.at( 1 ).dart() ).value() )
                };
                const std::array<size_t, 3> corner_ids{
                    vertex_ii( *mAtlas, vertices[0] ),
                    vertex_ii( *mAtlas, vertices[1] ),
                    vertex_ii( *mAtlas, vertices[2] )
                };

                const std::array<double, 3> bary{
                    expanded_coords( vertex_ii( *mAtlas, vertices[0] ) ),
                    expanded_coords( vertex_ii( *mAtlas, vertices[1] ) ),
                    expanded_coords( vertex_ii( *mAtlas, vertices[2] ) )
                };

                if( pt.mBaryCoordIsZero.at( corner_ids[1] ) )
                    return bary[0] * mTriMapping.vertPositions()( vertices[0] ) +
                           bary[2] * mTriMapping.vertPositions()( vertices[2] );

                const Eigen::Vector2d edge_midpoint =
                    0.5 * mTriMapping.vertPositions()( vertices[0] ) + 0.5 * mTriMapping.vertPositions()( vertices[2] );

                if( expanded_coords( corner_ids[ 0 ] ) > expanded_coords( corner_ids[ 2 ] ) )
                {
                    // right side triangle
                    const double theta0 = mBoundaryAngles.at( vert_ids( vertices[0] ) );
                    const double theta1 = mBoundaryAngles.at( vert_ids( vertices[1] ) );
                    return triangleWithArcEdgeEval( { bary[ 0 ] - bary[ 2 ], bary[ 1 ], 2 * bary[ 2 ] }, {theta0, theta1}, edge_midpoint );
                }
                else
                {
                    // left side triangle
                    const double theta1 = mBoundaryAngles.at( vert_ids( vertices[1] ) );
                    const double theta2 = mBoundaryAngles.at( vert_ids( vertices[2] ) );
                    return triangleWithArcEdgeEval( { bary[ 1 ], bary[ 2 ] - bary[ 0 ], 2 * bary[ 0 ] }, {theta1, theta2}, edge_midpoint );
                }
            }
            default:
                throw std::runtime_error( "Face with three or more boundary edges cannot exist in a TriangleMeshCircleMap");
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

    // See https://stackoverflow.com/a/23550032
    bool isBetweenAngles( double a, double b, const double test_angle )
    {
        a -= test_angle;
        b -= test_angle;
        a = util::normalizeAngle( a );
        b = util::normalizeAngle( b );
        if( a * b > 0 ) return false;
        return std::abs( a - b ) < std::numbers::pi;
    }

    std::optional<Eigen::Vector3d> invertTriangleWithArcEdgeEval( const Eigen::Vector2d& pt, const Eigen::Vector2d& non_bdry_point, const std::array<double, 2>& thetas )
    {
        const Eigen::Vector2d bdry_point = lineCircleIntersection( non_bdry_point, pt );
        const double bdry_point_theta = atan2( bdry_point( 1 ), bdry_point( 0 ) );

        if( not isBetweenAngles( thetas[ 0 ], thetas[ 1 ], bdry_point_theta ) ) return std::nullopt;

        const auto [ lambda_non_bdry, lambda_bdry ] = [&](){
            if( std::abs( bdry_point( 0 ) - non_bdry_point( 0 ) ) > std::abs( bdry_point( 1 ) - non_bdry_point( 1 ) ) )
                return inverseLinear( non_bdry_point(0), bdry_point(0), pt(0) );
            else
                return inverseLinear( non_bdry_point(1), bdry_point(1), pt(1) );
        }();

        const auto [lambda_v0, lambda_v1] = inverseLinear(
            util::normalizeAngle( thetas[ 0 ] - bdry_point_theta ), util::normalizeAngle( thetas[ 1 ] - bdry_point_theta ), 0.0 );

        Eigen::Vector3d bary = Eigen::Vector3d::Zero();
        bary( 0 ) = lambda_v0 * lambda_bdry;
        bary( 1 ) = lambda_v1 * lambda_bdry;
        bary( 2 ) = lambda_non_bdry;
        return bary;
    }

    std::optional<param::ParentPoint> TriangleMeshCircleMapping::maybeInverse( const topology::Face& f, const Eigen::Vector2d& pt ) const
    {
        const SmallVector<topology::Edge, 3> boundary_edges = maybeBoundaryEdges( mAtlas->cmap(), f );
        switch( boundary_edges.size() )
        {
            case 0: return mTriMapping.maybeInverse( f, pt );
            case 1:
            {
                const std::array<topology::Vertex, 3> vertices{
                    topology::Vertex( boundary_edges.at( 0 ).dart() ),
                    topology::Vertex( phi( mAtlas->cmap(), 1, boundary_edges.at( 0 ).dart() ).value() ),
                    topology::Vertex( phi( mAtlas->cmap(), -1, boundary_edges.at( 0 ).dart() ).value() )
                };
                
                const Eigen::Vector2d non_bdry_point = mTriMapping.vertPositions()( vertices[ 2 ] );

                const Eigen::Vector2d one_bdry_vert_point = mTriMapping.vertPositions()( vertices[ 0 ] );
                if( ( pt - non_bdry_point ).norm() < 1e-10 * ( one_bdry_vert_point - non_bdry_point ).norm() ) return mAtlas->parentPoint( vertices[ 2 ] );

                const auto vert_ids = indexingOrError( mAtlas->cmap(), 0 );
                const double theta0 = mBoundaryAngles.at( vert_ids( vertices[ 0 ] ) );
                const double theta1 = mBoundaryAngles.at( vert_ids( vertices[ 1 ] ) );

                return invertTriangleWithArcEdgeEval( pt, non_bdry_point, { theta0, theta1 } )
                    .transform( [&]( const Eigen::Vector3d& bary ) {
                        Eigen::Vector3d reordered_bary;
                        reordered_bary( vertex_ii( *mAtlas, vertices[0] ) ) = bary( 0 );
                        reordered_bary( vertex_ii( *mAtlas, vertices[1] ) ) = bary( 1 );
                        reordered_bary( vertex_ii( *mAtlas, vertices[2] ) ) = bary( 2 );
                        const param::ParentDomain pd = mAtlas->parentDomain( f );
                        return compressCoordinates( pd, reordered_bary, 1e-5 );
                    } );
            }
            case 2:
            {
                const std::array<topology::Vertex, 3> vertices{
                    topology::Vertex( boundary_edges.at( 0 ).dart() ),
                    topology::Vertex( boundary_edges.at( 1 ).dart() ),
                    topology::Vertex( phi( mAtlas->cmap(), 1, boundary_edges.at( 1 ).dart() ).value() )
                };

                const Eigen::Vector2d non_bdry_point =
                    0.5 * mTriMapping.vertPositions()( vertices[0] ) + 0.5 * mTriMapping.vertPositions()( vertices[2] );

                const Eigen::Vector2d one_bdry_vert_point = mTriMapping.vertPositions()( vertices[ 0 ] );
                if( ( pt - non_bdry_point ).norm() < 1e-10 * ( one_bdry_vert_point - non_bdry_point ).norm() )
                    return average( mAtlas->parentPoint( vertices[2] ), mAtlas->parentPoint( vertices[0] ) );

                const auto vert_ids = indexingOrError( mAtlas->cmap(), 0 );
                const double theta0 = mBoundaryAngles.at( vert_ids( vertices[ 0 ] ) );
                const double theta1 = mBoundaryAngles.at( vert_ids( vertices[ 1 ] ) );
                const double theta2 = mBoundaryAngles.at( vert_ids( vertices[ 2 ] ) );

                bool first_side = true;
                return invertTriangleWithArcEdgeEval( pt, non_bdry_point, { theta0, theta1 } )
                    .or_else( [&]() -> std::optional<Eigen::Vector3d> {
                        first_side = false;
                        return invertTriangleWithArcEdgeEval( pt, non_bdry_point, { theta1, theta2 } );
                    } )
                    .transform( [&]( const Eigen::Vector3d& bary ) {
                        Eigen::Vector3d reordered_bary;
                        if( first_side )
                        {
                            reordered_bary( vertex_ii( *mAtlas, vertices[0] ) ) = bary( 0 ) + 0.5 * bary( 2 );
                            reordered_bary( vertex_ii( *mAtlas, vertices[1] ) ) = bary( 1 );
                            reordered_bary( vertex_ii( *mAtlas, vertices[2] ) ) = 0.5 * bary( 2 );
                        }
                        else
                        {
                            reordered_bary( vertex_ii( *mAtlas, vertices[0] ) ) = 0.5 * bary( 2 );
                            reordered_bary( vertex_ii( *mAtlas, vertices[1] ) ) = bary( 0 );
                            reordered_bary( vertex_ii( *mAtlas, vertices[2] ) ) = bary( 1 ) + 0.5 * bary( 2 );
                        }
                        const param::ParentDomain pd = mAtlas->parentDomain( f );
                        return compressCoordinates( pd, reordered_bary, 1e-5 );
                    } );
            }
            default:
                throw std::runtime_error( "Face with three or more boundary edges cannot exist in a TriangleMeshCircleMap");
        }
    }

    std::optional<std::pair<topology::Cell, param::ParentPoint>> TriangleMeshCircleMapping::maybeInverse( const Vector3dMax& pt ) const
    {
        if( pt.size() != 2 ) throw std::invalid_argument( "TriangleMeshCircleMapping::maybeInverse requires a 2D point" );

        if( pt.norm() > 1.0 + 1e-15 ) return std::nullopt;
        std::optional<std::pair<topology::Cell, param::ParentPoint>> out;
        for( const auto& [f, bb] : mBoundingBoxes )
        {
            if( not bb.contains( pt ) ) continue;
            out = maybeInverse( f, pt ).transform( [&f]( const param::ParentPoint& ppt ) {
                return std::pair<topology::Cell, param::ParentPoint>{ f, ppt };
            } );
            if( out.has_value() ) return out;
        }
        return out;
    }

}