#include <Tracing.hpp>
#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include <set>
#include <iomanip>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <TetMeshCombinatorialMap.hpp>
#include <SweepInput.hpp>
#include <Logging.hpp>

constexpr bool LOG_TRACING = 0;

std::optional<Eigen::Vector3d> intersectionOf( const Ray<3>& ray,
                                               const Triangle<3>& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal )
{
    const Eigen::Vector3d normal = maybe_normal.value_or( triangleNormal( tri ) );
    const double ray_scaling = normal.dot( tri.v1 - ray.start_pos ) / normal.dot( ray.dir );
    LOG( LOG_TRACING ) << "| | Ray scaling: " << ray_scaling << std::endl;
    if( ray_scaling < 0 ) return std::nullopt;
    const Eigen::Vector3d intersection_point = ray.start_pos + ray_scaling * ray.dir;

    const bool inside = normal.dot( ( tri.v2 - tri.v1 ).cross( intersection_point - tri.v1 ) ) >= 0 and
                        normal.dot( ( tri.v3 - tri.v2 ).cross( intersection_point - tri.v2 ) ) >= 0 and
                        normal.dot( ( tri.v1 - tri.v3 ).cross( intersection_point - tri.v3 ) ) >= 0;

    return inside ? std::optional<Eigen::Vector3d>( intersection_point ) : std::nullopt;
}

std::optional<TracePoint> traceRayOnTet( const cgogn::CMap3& map,
                                         const cgogn::CMap3::Volume& v,
                                         const Ray<3>& ray,
                                         const std::vector<Normal>& normals )
{
    // The face on the input dart is the location that we start from.
    // Check all the other faces for intersection with the ray.
    LOG( LOG_TRACING ) << "Tracing on tet " << v << " from ray " << ray.start_pos.transpose() << " -> "
                       << ray.dir.transpose() << std::endl;
    bool found_it = false;
    std::pair<cgogn::CMap3::Face, Eigen::Vector3d> out;
    foreach_incident_edge( map, cgogn::CMap3::Face( v.dart_ ), [&]( cgogn::CMap3::Edge e ) {
        LOG( LOG_TRACING ) << "| Edge " << e << std::endl;
        const cgogn::CMap3::Face adj_face( phi2( map, e.dart_ ) );
        const Triangle<3> tri = triangleOfFace( map, adj_face );
        LOG( LOG_TRACING ) << "| | Triangle<3>: " << tri.v1.transpose() << " | " << tri.v2.transpose() << " | "
                           << tri.v3.transpose() << std::endl;

        const std::optional<Eigen::Vector3d> intersection =
            intersectionOf( ray, tri, normals.at( index_of( map, adj_face ) ).get( adj_face.dart_ ) );

        if( intersection.has_value() )
        {
            LOG( LOG_TRACING ) << "| | Found intersection: " << intersection.value().transpose() << std::endl;
            out.first = cgogn::CMap3::Face( phi3( map, adj_face.dart_ ) );
            out.second = intersection.value();
            found_it = true;
            return false;
        }
        return true;
    } );
    if( not found_it )
    {
        LOG( LOG_TRACING ) << "| | NO INTERSECTION\n";
        return std::nullopt;
    }
    return out;
}

std::optional<TracePoint2> traceRayOnTet( const topology::TetMeshCombinatorialMap& map,
                                          const topology::Volume& v,
                                          const Ray<3>& ray,
                                          const std::vector<Normal>& normals )
{
    // The face on the input dart is the location that we start from.
    // Check all the other faces for intersection with the ray.
    LOG( LOG_TRACING ) << "Tracing on tet " << v.dart().id() << " from ray " << ray.start_pos.transpose() << " -> "
                       << ray.dir.transpose() << std::endl;
    bool found_it = false;
    std::pair<topology::Face, Eigen::Vector3d> out;
    iterateAdjacentCells( map, topology::Face( v.dart() ), 1, [&]( const topology::Edge& e ) {
        LOG( LOG_TRACING ) << "| Edge " << e.dart().id() << std::endl;
        const topology::Face adj_face( phi( map, 2, e.dart() ).value() );
        const Triangle<3> tri = triangleOfFace( map, adj_face );
        LOG( LOG_TRACING ) << "| | Triangle<3>: " << tri.v1.transpose() << " | " << tri.v2.transpose() << " | "
                           << tri.v3.transpose() << std::endl;

        const std::optional<Eigen::Vector3d> intersection =
            intersectionOf( ray, tri, normals.at( map.faceId( adj_face ) ).get( adj_face.dart() ) );

        if( intersection.has_value() )
        {
            LOG( LOG_TRACING ) << "| | Found intersection: " << intersection.value().transpose() << std::endl;
            const auto maybe_phi3 = phi( map, 3, adj_face.dart() );
            out.first = topology::Face( maybe_phi3.has_value() ? maybe_phi3.value() : adj_face.dart() ); // FIXME?
            out.second = intersection.value();
            found_it = true;
            return false;
        }
        return true;
    } );
    if( not found_it )
    {
        LOG( LOG_TRACING ) << "| | NO INTERSECTION\n";
        return std::nullopt;
    }
    return out;
}

void tracingDebugOutput( const cgogn::CMap3& map,
                         const cgogn::CMap3::Volume& v,
                         const Ray<3>& ray,
                         const SimplicialComplex& line,
                         SimplicialComplex& tets,
                         const size_t n )
{
    // Add the triangles of the current tet to tets
    foreach_incident_face( map, v, [&]( cgogn::CMap3::Face f ) {
        const Triangle<3> from_face = triangleOfFace( map, f );
        const size_t next_vid = tets.points.size();
        tets.simplices.push_back( { next_vid + 0, next_vid + 1, next_vid + 2 } );
        tets.points.push_back( from_face.v1 );
        tets.points.push_back( from_face.v2 );
        tets.points.push_back( from_face.v3 );
        return true;
    } );

    // Create a simplicial complex with just one face: the triangle it's coming from
    const Triangle<3> from_face = triangleOfFace( map, cgogn::CMap3::Face( v.dart_ ) );
    const SimplicialComplex from_face_complex( { { { 0, 1, 2 } }, { from_face.v1, from_face.v2, from_face.v3 } } );
    const io::VTKOutputObject from_face_output( from_face_complex );

    // Create a simplicial complex for the ray
    const SimplicialComplex ray_complex( { { { 0 } }, { ray.start_pos } } );
    io::VTKOutputObject ray_output( ray_complex );
    ray_output.addVertexField( "ray", ray.dir.transpose() );

    const io::VTKOutputObject line_output( line );
    const io::VTKOutputObject tets_output( tets );

    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << n;
    std::string n_str( ss.str() );

    io::outputSimplicialFieldToVTK( from_face_output, "from_face_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( ray_output, "ray_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( line_output, "line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( tets_output, "tets_" + n_str + ".vtu" );
}

void tracingDebugOutput( const topology::TetMeshCombinatorialMap& map,
                         const topology::Volume& v,
                         const Ray<3>& ray,
                         const SimplicialComplex& line,
                         SimplicialComplex& tets,
                         const size_t n )
{
    // Add the triangles of the current tet to tets
    iterateAdjacentCells( map, v, 2, [&]( const topology::Face& f ) {
        const Triangle<3> from_face = triangleOfFace( map, f );
        const size_t next_vid = tets.points.size();
        tets.simplices.push_back( { next_vid + 0, next_vid + 1, next_vid + 2 } );
        tets.points.push_back( from_face.v1 );
        tets.points.push_back( from_face.v2 );
        tets.points.push_back( from_face.v3 );
        return true;
    } );

    // Create a simplicial complex with just one face: the triangle it's coming from
    const Triangle<3> from_face = triangleOfFace( map, topology::Face( v.dart() ) );
    const SimplicialComplex from_face_complex( { { { 0, 1, 2 } }, { from_face.v1, from_face.v2, from_face.v3 } } );
    const io::VTKOutputObject from_face_output( from_face_complex );

    // Create a simplicial complex for the ray
    const SimplicialComplex ray_complex( { { { 0 } }, { ray.start_pos } } );
    io::VTKOutputObject ray_output( ray_complex );
    ray_output.addVertexField( "ray", ray.dir.transpose() );

    const io::VTKOutputObject line_output( line );
    const io::VTKOutputObject tets_output( tets );

    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << n;
    std::string n_str( ss.str() );

    io::outputSimplicialFieldToVTK( from_face_output, "from_face_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( ray_output, "ray_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( line_output, "line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( tets_output, "tets_" + n_str + ".vtu" );
}

void tracingDebugOutput( const cgogn::CMap3& map,
                         const cgogn::CMap3::Face& curr_face,
                         const double start_pos,
                         const Eigen::Ref<const Eigen::VectorXd> field,
                         const SimplicialComplex& line,
                         SimplicialComplex& tris,
                         std::vector<double>& vertex_values,
                         const size_t n )
{
    const Triangle<3> tri3d = triangleOfFace( map, curr_face );
    const auto& d = curr_face.dart_;
    const auto field_values = field( { index_of( map, cgogn::CMap3::Vertex( d ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) ) } );

    // Add the triangles of the current tet to tets
    const size_t next_vid = tris.points.size();
    tris.simplices.push_back( { next_vid + 0, next_vid + 1, next_vid + 2 } );
    tris.points.push_back( tri3d.v1 );
    tris.points.push_back( tri3d.v2 );
    tris.points.push_back( tri3d.v3 );
    vertex_values.push_back( field_values( 0 ) );
    vertex_values.push_back( field_values( 1 ) );
    vertex_values.push_back( field_values( 2 ) );

    // Create a simplicial complex with just one face: the triangle it's coming from
    const SimplicialComplex from_line_complex( { { { 0, 1 } }, { tri3d.v1, tri3d.v2 } } );
    const io::VTKOutputObject from_line_output( from_line_complex );

    // Create a simplicial complex for the ray
    const Ray<3> ray( { ( 1.0 - start_pos ) * tri3d.v1 + start_pos * tri3d.v2, gradient( tri3d, field_values ) } );
    const SimplicialComplex ray_complex( { { { 0 } }, { ray.start_pos } } );
    io::VTKOutputObject ray_output( ray_complex );
    ray_output.addVertexField( "ray", ray.dir.transpose() );

    const io::VTKOutputObject line_output( line );
    io::VTKOutputObject tets_output( tris );
    tets_output.addVertexField( "field", Eigen::Map<Eigen::VectorXd>( vertex_values.data(), vertex_values.size() ) );

    std::stringstream ss;
    ss << std::setw( 3 ) << std::setfill( '0' ) << n;
    std::string n_str( ss.str() );

    io::outputSimplicialFieldToVTK( from_line_output, "bdry_from_line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( ray_output, "bdry_ray_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( line_output, "bdry_line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( tets_output, "bdry_tris_" + n_str + ".vtu" );
}

std::optional<TracePoint> traceFaceAverageField( const cgogn::CMap3& map,
                                                 const cgogn::CMap3::Face& f,
                                                 const Eigen::Vector3d& start_point,
                                                 const Eigen::MatrixX3d& field,
                                                 const std::vector<Normal>& normals )
{
    const cgogn::CMap3::Face f_opp( phi3( map, f.dart_ ) );
    const cgogn::CMap3::Volume vol( f.dart_ );
    const cgogn::CMap3::Volume vol_opp( f_opp.dart_ );
    const auto average_field = 0.5 * ( field.row( index_of( map, vol ) ) + field.row( index_of( map, vol_opp ) ) );
    const double field_dir = normals.at( index_of( map, f ) ).get( f.dart_ ).dot( average_field );
    if( field_dir > 0 )
        return traceRayOnTet( map, vol, { start_point, average_field }, normals );
    else
        return traceRayOnTet( map, vol_opp, { start_point, average_field }, normals );
}

std::optional<TracePoint2> traceFaceAverageField( const topology::TetMeshCombinatorialMap& map,
                                                  const topology::Face& f,
                                                  const Eigen::Vector3d& start_point,
                                                  const Eigen::MatrixX3d& field,
                                                  const std::vector<Normal>& normals )
{
    LOG( LOG_TRACING ) << "| | | Tracing on face average field\n";
    return phi( map, 3, f.dart() ).and_then( [&]( const topology::Dart& d ) {
        LOG( LOG_TRACING ) << "| | | | Has phi3\n";
        const topology::Face f_opp( d );
        const topology::Volume vol( f.dart() );
        const topology::Volume vol_opp( d );
        const auto average_field = 0.5 * ( field.row( map.elementId( vol ) ) + field.row( map.elementId( vol_opp ) ) );
        const double field_dir = normals.at( map.faceId( f ) ).get( f.dart() ).dot( average_field );
        if( field_dir > 0 )
            return traceRayOnTet( map, vol, { start_point, average_field }, normals );
        else
            return traceRayOnTet( map, vol_opp, { start_point, average_field }, normals );
    } );
}

SimplicialComplex traceField( const cgogn::CMap3& map,
                              const cgogn::CMap3::Face& f,
                              const Eigen::Vector3d& start_point,
                              const Eigen::MatrixX3d& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output )
{
    cgogn::CMap3::Face curr_face = f;
    Eigen::Vector3d curr_point = start_point;

    SimplicialComplex complex;
    SimplicialComplex debug_tets;
    size_t n = 0;

    complex.points.push_back( curr_point );
    while( not is_boundary( map, curr_face.dart_ ) )
    {
        const cgogn::CMap3::Volume curr_vol( curr_face.dart_ );
        const Ray<3> search_ray( { curr_point, field.row( index_of( map, curr_vol ) ) } );
        if( debug_output ) tracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n++ );
        std::optional<TracePoint> next_point = traceRayOnTet( map, curr_vol, search_ray, normals );
        if( not next_point.has_value() )
        {
            next_point = traceFaceAverageField( map, cgogn::CMap3::Face( curr_vol.dart_ ), curr_point, field, normals );
            if( not next_point.has_value() ) throw( std::runtime_error( "Untraceable field" ) );
        }
        curr_face = next_point.value().first;
        curr_point = next_point.value().second;
        complex.points.push_back( curr_point );
        complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
    }

    return complex;
}

SimplicialComplex traceField( const topology::TetMeshCombinatorialMap& map,
                              const topology::Face& f,
                              const Eigen::Vector3d& start_point,
                              const Eigen::MatrixX3d& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output )
{
    topology::Face curr_face = f;
    Eigen::Vector3d curr_point = start_point;

    SimplicialComplex complex;
    SimplicialComplex debug_tets;
    size_t n = 0;

    complex.points.push_back( curr_point );
    LOG( LOG_TRACING ) << "Starting a trace\n";
    do
    {
        const topology::Volume curr_vol( curr_face.dart() );
        const Ray<3> search_ray( { curr_point, field.row( map.elementId( curr_vol ) ) } );
        if( debug_output ) tracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n++ );
        std::optional<TracePoint2> next_point = traceRayOnTet( map, curr_vol, search_ray, normals );
        if( not next_point.has_value() )
        {
            next_point = traceFaceAverageField( map, topology::Face( curr_vol.dart() ), curr_point, field, normals );
            if( not next_point.has_value() )
            {
                throw( std::runtime_error( "Untraceable field" ) );
            }
        }
        curr_face = next_point.value().first;
        curr_point = next_point.value().second;
        complex.points.push_back( curr_point );
        complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
    } while( not onBoundary( map, curr_face.dart() ) );

    return complex;
}

std::optional<double> barycentricIntersectionOf( const Ray<2>& ray, const Segment<2>& line )
{
    const Ray<2> line_as_ray( { line.start_pos, line.end_pos - line.start_pos } );
    const Eigen::Ref<const Eigen::Vector2d> start_pos_diff = line.start_pos - ray.start_pos;

    const double ray_multiplier_numerator =
        start_pos_diff( 1 ) * line_as_ray.dir( 0 ) - start_pos_diff( 0 ) * line_as_ray.dir( 1 );

    const double line_multiplier_numerator = start_pos_diff( 1 ) * ray.dir( 0 ) - start_pos_diff( 0 ) * ray.dir( 1 );

    const double line_multiplier_denominator =
        line_as_ray.dir( 0 ) * ray.dir( 1 ) - line_as_ray.dir( 1 ) * ray.dir( 0 );

    if( ray_multiplier_numerator * line_multiplier_denominator < 0 ) return std::nullopt;
    if( line_multiplier_numerator * line_multiplier_denominator < 0 ) return std::nullopt;
    if( std::abs( line_multiplier_denominator ) < std::abs( line_multiplier_numerator ) ) return std::nullopt;

    const double u = line_multiplier_numerator / line_multiplier_denominator;

    return u;
}

std::optional<Eigen::Vector2d> intersectionOf( const Ray<2>& ray, const Segment<2>& line )
{
    return barycentricIntersectionOf( ray, line ).and_then( [&]( const double& u ) -> std::optional<Eigen::Vector2d> {
        return ( 1 - u ) * line.start_pos + u * line.end_pos;
    } );
}

constexpr bool FORWARD = true;
constexpr bool BACKWARD = false;
std::optional<std::pair<bool, double>> traceGradientOnTri( const Triangle<3>& tri3d,
                                                           const double edge_barycentric_coord,
                                                           const Eigen::Ref<const Eigen::Vector3d> field_values )
{
    // move the triangle into the xy plane
    const Eigen::Vector3d e1 = ( tri3d.v2 - tri3d.v1 ).normalized();
    const Eigen::Vector3d e2 = ( tri3d.v3 - tri3d.v1 - e1.dot( tri3d.v3 - tri3d.v1 ) * e1 ).normalized();

    const Eigen::Vector2d v1_2d( 0, 0 );
    const Eigen::Vector2d v2_2d( ( tri3d.v2 - tri3d.v1 ).norm(), 0 );
    const Eigen::Vector2d v3_2d( e1.dot( tri3d.v3 - tri3d.v1 ), e2.dot( tri3d.v3 - tri3d.v1 ) );

    // calculate the gradient in the xy plane
    const double& f1 = field_values( 0 );
    const double& f2 = field_values( 1 );
    const double& f3 = field_values( 2 );

    const Eigen::Rotation2Dd rot90( std::numbers::pi / 2.0 );
    const double twice_area = v2_2d( 0 ) * v3_2d( 1 );

    const auto grad_s_i = [&]( const Segment<2>& edge_i ) -> Eigen::Vector2d {
        const auto edge_diff = edge_i.end_pos - edge_i.start_pos;
        return 1.0 / twice_area * ( rot90 * edge_diff );
    };

    const Eigen::Vector2d gradient =
        f1 * grad_s_i( { v2_2d, v3_2d } ) + f2 * grad_s_i( { v3_2d, v1_2d } ) + f3 * grad_s_i( { v1_2d, v2_2d } );

    // iterate the edges and perform line ray intersections
    const Ray<2> ray( { edge_barycentric_coord * v2_2d, gradient } );

    const auto run_intersection = [&]( const bool i, const Segment<2>& line ) {
        return barycentricIntersectionOf( ray, line )
            .and_then( [&]( const double& u ) -> std::optional<std::pair<unsigned int, double>> {
                return std::pair<bool, double>{ i, u };
            } );
    };

    return run_intersection( FORWARD, { v2_2d, v3_2d } ).or_else( [&]() {
        return run_intersection( BACKWARD, { v3_2d, v1_2d } );
    } );
}

std::optional<std::pair<cgogn::CMap3::Edge, double>> traceGradientOnTri( const cgogn::CMap3& map,
                                                                         const cgogn::CMap3::Face& f,
                                                                         const double edge_barycentric_coord,
                                                                         const Eigen::VectorXd& field_values )
{
    const Triangle<3> tri3d = triangleOfFace( map, f );

    // calculate the gradient in the xy plane
    const cgogn::Dart& d = f.dart_;
    const auto field = field_values( { index_of( map, cgogn::CMap3::Vertex( d ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) ) } );

    return traceGradientOnTri( tri3d, edge_barycentric_coord, field )
        .and_then( [&]( const std::pair<bool, double>& pr ) {
            const cgogn::CMap3::Edge e( pr.first ? phi<1, 2>( map, f.dart_ ) : phi<-1, 2>( map, f.dart_ ) );
            return std::optional<std::pair<cgogn::CMap3::Edge, double>>( { e, 1.0 - pr.second } );
        } );
}

std::optional<std::pair<cgogn::CMap3::Edge, double>> traceGradientOnTriPair( const cgogn::CMap3& map,
                                                                             const cgogn::CMap3::Face& f,
                                                                             const double edge_barycentric_coord,
                                                                             const Eigen::VectorXd& field_values )
{
    const Triangle<3> tri3d = triangleOfFace( map, f );
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const Eigen::Vector3d& opp_v =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( phi<2, -1>( map, f.dart_ ) ) );

    // calculate the gradient in the xy plane
    const cgogn::Dart& d = f.dart_;
    const auto field = field_values( { index_of( map, cgogn::CMap3::Vertex( d ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) ),
                                       index_of( map, cgogn::CMap3::Vertex( cgogn::phi<2, -1>( map, d ) ) ) } );

    // move the triangle into the xy plane
    const Eigen::Vector3d e1 = ( tri3d.v2 - tri3d.v1 ).normalized();
    const Eigen::Vector3d e2 = ( tri3d.v3 - tri3d.v1 - e1.dot( tri3d.v3 - tri3d.v1 ) * e1 ).normalized();
    const Eigen::Vector3d e2_prime = -( opp_v - tri3d.v1 - e1.dot( opp_v - tri3d.v1 ) * e1 ).normalized();

    const Eigen::Vector2d v1_2d( 0, 0 );
    const Eigen::Vector2d v2_2d( ( tri3d.v2 - tri3d.v1 ).norm(), 0 );
    const Eigen::Vector2d v3_2d( e1.dot( tri3d.v3 - tri3d.v1 ), e2.dot( tri3d.v3 - tri3d.v1 ) );
    const Eigen::Vector2d v4_2d( e1.dot( opp_v - tri3d.v1 ), e2_prime.dot( opp_v - tri3d.v1 ) );

    // calculate the gradient in the xy plane
    const double& f1 = field( 0 );
    const double& f2 = field( 1 );
    const double& f3 = field( 2 );
    const double& f4 = field( 3 );

    const Eigen::Rotation2Dd rot90( std::numbers::pi / 2.0 );
    const double twice_area_1 = v2_2d( 0 ) * v3_2d( 1 );
    const double twice_area_2 = -v2_2d( 0 ) * v4_2d( 1 );

    const auto grad_s_i = [&rot90]( const Segment<2>& edge_i, const double& twice_area ) -> Eigen::Vector2d {
        const auto edge_diff = edge_i.end_pos - edge_i.start_pos;
        return 1.0 / twice_area * ( rot90 * edge_diff );
    };

    const Eigen::Vector2d gradient_1 = f1 * grad_s_i( { v2_2d, v3_2d }, twice_area_1 ) +
                                       f2 * grad_s_i( { v3_2d, v1_2d }, twice_area_1 ) +
                                       f3 * grad_s_i( { v1_2d, v2_2d }, twice_area_1 );

    const Eigen::Vector2d gradient_2 = f1 * grad_s_i( { v4_2d, v2_2d }, twice_area_2 ) +
                                       f2 * grad_s_i( { v1_2d, v4_2d }, twice_area_2 ) +
                                       f4 * grad_s_i( { v2_2d, v1_2d }, twice_area_2 );

    const Eigen::Vector2d gradient = ( gradient_1 + gradient_2 ) * 0.5;

    // iterate the edges and perform line ray intersections
    const Ray<2> ray( { edge_barycentric_coord * v2_2d, gradient } );

    const auto run_intersection = [&]( const bool i, const Segment<2>& line ) {
        return barycentricIntersectionOf( ray, line )
            .and_then( [&]( const double& u ) -> std::optional<std::pair<unsigned int, double>> {
                return std::pair<bool, double>{ i, u };
            } );
    };

    if( gradient( 1 ) >= 0 )
    {
        return run_intersection( FORWARD, { v2_2d, v3_2d } )
            .or_else( [&]() {
                return run_intersection( BACKWARD, { v3_2d, v1_2d } );
            } )
            .and_then( [&]( const std::pair<bool, double>& pr ) {
                const cgogn::CMap3::Edge e( pr.first ? phi<1, 2>( map, f.dart_ ) : phi<-1, 2>( map, f.dart_ ) );
                return std::optional<std::pair<cgogn::CMap3::Edge, double>>( { e, 1.0 - pr.second } );
            } );
    }
    else
    {
        return run_intersection( FORWARD, { v4_2d, v2_2d } )
            .or_else( [&]() {
                return run_intersection( BACKWARD, { v1_2d, v4_2d } );
            } )
            .and_then( [&]( const std::pair<bool, double>& pr ) {
                const cgogn::CMap3::Edge e( pr.first ? phi<2, -1, 2>( map, f.dart_ ) : phi<2, 1, 2>( map, f.dart_ ) );
                return std::optional<std::pair<cgogn::CMap3::Edge, double>>( { e, 1.0 - pr.second } );
            } );
    }
}

SimplicialComplex traceBoundaryField( const cgogn::CMap3& map,
                                      const cgogn::CMap3::Edge& e,
                                      const double& start_point,
                                      const Eigen::VectorXd& field,
                                      const std::vector<bool>& target_vertex_marker,
                                      const bool debug_output,
                                      const std::function<void( const cgogn::CMap3::Face& )>& face_callback )
{
    const auto is_target = [&map, &target_vertex_marker]( const cgogn::CMap3::Edge& e ) {
        return target_vertex_marker.at( index_of( map, cgogn::CMap3::Vertex( e.dart_ ) ) ) and
               target_vertex_marker.at( index_of( map, cgogn::CMap3::Vertex( phi1( map, e.dart_ ) ) ) );
    };

    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const auto expand_barycentric = [&map, &position]( const cgogn::CMap3::Edge& e, const double coord ) {
        const cgogn::Dart& d = e.dart_;
        const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
        const Eigen::Vector3d& pos2 =
            cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );

        return ( 1.0 - coord ) * pos1 + coord * pos2;
    };

    cgogn::CMap3::Edge curr_edge = e;
    double curr_point = start_point;

    SimplicialComplex traced_line;
    SimplicialComplex debug_tris;
    std::vector<double> debug_field_values;
    size_t n = 0;

    traced_line.points.push_back( expand_barycentric( curr_edge, curr_point ) );
    while( not is_target( curr_edge ) )
    {
        const cgogn::CMap3::Face curr_face( curr_edge.dart_ );
        face_callback( curr_face );
        if( debug_output )
            tracingDebugOutput( map, curr_face, curr_point, field, traced_line, debug_tris, debug_field_values, n++ );
        std::optional<std::pair<cgogn::CMap3::Edge, double>> next_point =
            traceGradientOnTri( map, curr_face, curr_point, field );

        if( not next_point.has_value() )
        {
            next_point = traceGradientOnTriPair( map, curr_face, curr_point, field );
            if( not next_point.has_value() )
            {
                throw( "No intersection on boundary tracing" );
                break;
            }
        }

        curr_edge = next_point.value().first;
        curr_point = next_point.value().second;
        traced_line.points.push_back( expand_barycentric( curr_edge, curr_point ) );
        traced_line.simplices.push_back( { traced_line.points.size() - 2, traced_line.points.size() - 1 } );
    }

    return traced_line;
}