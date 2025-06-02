#include <OrbifoldMapping.hpp>

using namespace mapping;

enum class RotationType
{
    None,
    Clockwise,
    CounterClockwise,
    OneEighty
};

Eigen::Vector2d transform( const Eigen::Vector2d& pt, const Eigen::Vector2i& translation, const RotationType& rotation )
{
    Eigen::Vector2d out = pt - translation.cast<double>();
    switch( rotation )
    {
        case RotationType::None:
            break;
        case RotationType::Clockwise:
            out = Eigen::Vector2d( out( 1 ), 1 - out( 0 ) );
            break;
        case RotationType::CounterClockwise:
            out = Eigen::Vector2d( 1 - out( 1 ), out( 0 ) );
            break;
        case RotationType::OneEighty:
            out = Eigen::Vector2d( 1 - out( 0 ), 1 - out( 1 ) );
            break;
    }
    return out;
}

std::pair<Eigen::Vector2i, RotationType> getTransform( const Eigen::Vector2d& pt )
{
    const Eigen::Vector2i int_part = pt.array().floor().cast<int>();

    if( int_part( 0 ) % 2 != 0 and int_part( 1 ) % 2 != 0 )
        return { int_part, RotationType::OneEighty };
    else if( int_part( 1 ) % 2 != 0 )
        return { int_part, RotationType::CounterClockwise };
    else if( int_part( 0 ) % 2 != 0 )
        return { int_part, RotationType::Clockwise };
    else
        return { int_part, RotationType::None };
}

Eigen::Vector2d transformToCanonicalDomain( const Eigen::Vector2d& pt )
{
    const auto [translation, rotation] = getTransform( pt );
    return transform( pt, translation, rotation );
}

OrbifoldMapping::OrbifoldMapping( const std::shared_ptr<const param::TriangleParametricAtlas>& atlas,
                                  const VertexPositionsFunc& vertex_positions )
    : mTriangleMapping( atlas, vertex_positions, 2, false )
{
    iterateCellsWhile( atlas->cmap(), 2, [&]( const topology::Face& f ) {
        const topology::Vertex v = originVertex( *atlas, f );
        const Triangle<2> tri = triangleOfFace<2>( atlas->cmap(), vertex_positions, topology::Face( v.dart() ) );
        const AABB aabb = aabbFromTriangle( tri );

        // If the triangle max is greater than 1 in either dimension, or min is less than 0 in either dimension,
        // we need to create copies of it. If the whole triangle is outside the canonical domain, we can just use
        // the rotated copy.
        const auto add_with_rotation = [&]( const Eigen::Vector2i& translation,
                                            const RotationType& rotation ) {
            const Triangle<2> tri_copy{ transform( tri.v1, translation, rotation ),
                                        transform( tri.v2, translation, rotation ),
                                        transform( tri.v3, translation, rotation ) };
            const AABB aabb_copy = aabbFromTriangle( tri_copy );
            mTriangles.emplace_back( OrbifoldTriangle{
                f,
                { std::pair<AABB, Triangle<2>>( aabb, tri ), std::pair<AABB, Triangle<2>>( aabb_copy, tri_copy ) } } );
        };
        if( ( aabb.min().array() > 1.0 ).any() or ( aabb.max().array() < 0.0 ).any() )
        {
            // Do we need to check the 180 degrees thing here as well?
            const Triangle<2> tri_copy{ transformToCanonicalDomain( tri.v1 ),
                                        transformToCanonicalDomain( tri.v2 ),
                                        transformToCanonicalDomain( tri.v3 ) };
            const AABB aabb_copy = aabbFromTriangle( tri_copy );
            mTriangles.push_back( OrbifoldTriangle{ f, { std::pair<AABB, Triangle<2>>( aabb_copy, tri_copy ) } } );
        }
        else if( ( aabb.max().array() > 1.0 ).any() )
        {
            const auto [translation, rotation] = [&]() -> std::pair<Eigen::Vector2i, RotationType> {
                const auto [translation, rotation] = getTransform( aabb.max() );
                if( rotation == RotationType::OneEighty )
                {
                    Eigen::Index max_index;
                    const double max_val = aabb.max().maxCoeff( &max_index );
                    const Eigen::Vector2d point = tri.v1( max_index ) == max_val ? tri.v1 : ( tri.v2( max_index ) == max_val ? tri.v2 : tri.v3 );
                    return getTransform( point );
                }
                else return { translation, rotation };
            }();
            add_with_rotation( translation, rotation );
        }
        else if( ( aabb.min().array() < 0.0 ).any() )
        {
            // Do we need to check the 180 degrees thing here as well?
            const auto [translation, rotation] = getTransform( aabb.min() );
            add_with_rotation( translation, rotation );
        }
        else
        {
            // The triangle is in the canonical domain, so we can just use it
            mTriangles.push_back( OrbifoldTriangle{ f, { std::pair<AABB, Triangle<2>>( aabb, tri ) } } );
        }

        return true;
    } );
}

const param::TriangleParametricAtlas& OrbifoldMapping::parametricAtlas() const
{
    return mTriangleMapping.parametricAtlas();
}

Eigen::VectorXd OrbifoldMapping::evaluate( const topology::Cell& c, const param::ParentPoint& pt ) const
{
    return transformToCanonicalDomain( mTriangleMapping.evaluate( c, pt ).head( 2 ) );
}

size_t OrbifoldMapping::spatialDim() const
{
    return 2;
}

std::optional<std::pair<topology::Cell, param::ParentPoint>> OrbifoldMapping::maybeInverse( const Vector3dMax& pt ) const
{
    const Eigen::Vector2d canonical_pt = transformToCanonicalDomain( pt.head( 2 ) );
    for( const auto& tri : mTriangles )
    {
        for( const auto& [bb, t] : tri.triangles )
        {
            if( bb.contains( canonical_pt ) )
            {
                const param::ParentDomain pd = parametricAtlas().parentDomain( tri.face );
                const auto maybe_bary = invertTriangleMap( t, canonical_pt );
                if( maybe_bary.has_value() )
                {
                    const param::ParentPoint ppt = compressCoordinates( pd, maybe_bary.value(), 1e-5 );
                    return std::pair<topology::Cell, param::ParentPoint>{ tri.face, ppt };
                }
            }
        }
    }
    throw std::runtime_error( "No inverse found on plane-tiling orbifold " + std::to_string( pt( 0 ) ) + " " +
                              std::to_string( pt( 1 ) ) );
}

std::pair<topology::Cell, param::ParentPoint> OrbifoldMapping::closestPoint( const Vector3dMax& pt ) const
{
    return maybeInverse( pt ).value();
}

const VertexPositionsFunc& OrbifoldMapping::vertPositions() const
{
    return mTriangleMapping.vertPositions();
}
