#include <Tracing.hpp>
#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include <set>
#include <iomanip>
#include <TetMeshCombinatorialMap.hpp>
#include <SweepInput.hpp>
#include <Logging.hpp>
#include <GlobalCellMarker.hpp>

namespace reparam
{

constexpr bool LOG_TRACING = 0;

////////////////////////////////////////
////////////////////////////////////////
/// GENERAL UTILITIES
////////////////////////////////////////

static constexpr bool LOG_TRACING_CELL = false;
std::optional<topology::Cell> tracingStartCell( const topology::CombinatorialMap& map,
                                                const topology::Cell& possible_cell,
                                                const std::function<Eigen::Vector3d( const topology::Cell& )>& normal,
                                                const std::function<Eigen::Vector3d( const topology::Cell& )>& grad )
{
    LOG( LOG_TRACING_CELL ) << "Searching for tracing cell...\n";
    const size_t dim = map.dim();

    LOG( LOG_TRACING_CELL ) << "Cell dim: " << possible_cell.dim() << std::endl;
    std::optional<topology::Cell> out;
    iterateAdjacentCells( map, possible_cell, dim, [&]( const topology::Cell& elem ) {
        const Eigen::Vector3d this_grad = grad( elem );
        LOG( LOG_TRACING_CELL ) << "Grad: " << this_grad.transpose() << std::endl;
        const bool found_cell =
            iterateAdjacentCellsOfRestrictedCell( map,
                                                  topology::Cell( elem.dart(), possible_cell.dim() ),
                                                  elem.dim(),
                                                  dim - 1,
                                                  [&]( const topology::Cell& interface ) {
                                                      const Eigen::Vector3d this_normal = normal( interface );
                                                      LOG( LOG_TRACING_CELL ) << "Normal: " << this_normal.transpose() << std::endl;
                                                      LOG( LOG_TRACING_CELL ) << "Dot product: " << this_normal.dot( this_grad ) << std::endl;
                                                      if( this_normal.dot( this_grad ) < 0 ) return false;
                                                      return true;
                                                  } );
        if( found_cell )
        {
            out.emplace( topology::Cell( elem.dart(), possible_cell.dim() ) );
            return false;
        }
        return true;
    } );
    return out;
}

std::optional<topology::Cell>
    tracingAverageStartCell( const topology::CombinatorialMap& map,
                             const topology::Cell& possible_cell,
                             const std::function<Eigen::Vector3d( const topology::Cell& )>& normal,
                             const std::function<Eigen::Vector3d( const topology::Cell& )>& grad,
                             const std::function<Eigen::Vector3d( const topology::Cell& )>& edge_func )
{
    LOG( LOG_TRACING_CELL ) << "Searching for tracing cell...\n";
    const size_t dim = map.dim();

    const auto vertex_ids = indexingOrError( map, 0 );

    LOG( LOG_TRACING_CELL ) << "Cell dim: " << possible_cell.dim() << std::endl;
    std::optional<topology::Cell> out;
    iterateAdjacentCells( map, possible_cell, dim - 1, [&]( const topology::Cell& I ) {
        const auto d_opp = phi( map, dim, I.dart() );
        if( not d_opp.has_value() ) return true;

        const Eigen::Vector3d grad_left = grad( topology::Cell( I.dart(), dim ) );
        const Eigen::Vector3d grad_right = grad( topology::Cell( d_opp.value(), dim ) );

        const Eigen::Vector3d normal_left = normal( I );
        const Eigen::Vector3d normal_right = normal( topology::Cell( d_opp.value(), I.dim() ) );

        const Eigen::Vector3d edge = edge_func( I );

        if( ( normal_left.dot( grad_left ) > 0 ) == ( normal_right.dot( grad_right ) > 0 ) )
        {
            const bool aligned = vertex_ids( topology::Vertex( I.dart() ) ) == vertex_ids( possible_cell );
            const bool grads_along_edge = ( edge.dot( grad_left ) + edge.dot( grad_right ) ) >= 0 == aligned;
            if( grads_along_edge )
            {
                out.emplace( topology::Cell( aligned ? I.dart() : d_opp.value(), possible_cell.dim() ) );
                return false;
            }
        }
        return true;
    } );
    return out;
}

////////////////////////////////////////
////////////////////////////////////////
/// DEBUG
////////////////////////////////////////

void interiorTracingDebugOutput( const topology::TetMeshCombinatorialMap& map,
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

void boundaryTracingDebugOutput( const topology::CombinatorialMap& map,
                                 const VertexPositionsFunc& positions,
                                 const topology::Cell& curr_cell,
                                 const double start_pos,
                                 const Eigen::Ref<const Eigen::VectorXd> field,
                                 const SimplicialComplex& line,
                                 SimplicialComplex& tris,
                                 std::vector<double>& vertex_values,
                                 const size_t n )
{
    const Triangle<3> tri3d = triangleOfFace<3>( map, positions, topology::Face( curr_cell.dart() ) );

    const auto& d = curr_cell.dart();
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto field_values = field( { vertex_ids( topology::Vertex( d ) ),
                                       vertex_ids( topology::Vertex( phi( map, 1, d ).value() ) ),
                                       vertex_ids( topology::Vertex( phi( map, -1, d ).value() ) ) } );

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
    const SimplicialComplex from_line_complex = [&](){
        if( curr_cell.dim() == 1 ) return SimplicialComplex( { { { 0, 1 } }, { tri3d.v1, tri3d.v2 } } );
        else return SimplicialComplex( { { { 0 } }, { tri3d.v1 } } );
    }();
    const io::VTKOutputObject from_line_output( from_line_complex );

    // Create a simplicial complex for the ray
    const auto grad = gradient( tri3d, field_values );
    const Eigen::Vector3d curr_pos = [&]() -> Eigen::Vector3d {
        if( curr_cell.dim() == 1 ) return ( 1.0 - start_pos ) * tri3d.v1 + start_pos * tri3d.v2;
        else return tri3d.v1;
    }();
    const Ray<3> ray( { curr_pos, grad } );
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

////////////////////////////////////////
////////////////////////////////////////
/// INTERIOR TRACING
////////////////////////////////////////

std::optional<Eigen::Vector3d> intersectionOf( const Ray<3>& ray,
                                               const Triangle<3>& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal )
{
    const Eigen::Vector3d normal = maybe_normal.has_value() ? maybe_normal.value() : triangleNormal( tri );
    const double ray_scaling = normal.dot( tri.v1 - ray.start_pos ) / normal.dot( ray.dir );
    LOG( LOG_TRACING ) << "| | Ray scaling: " << ray_scaling << std::endl;
    if( ray_scaling < 0 ) return std::nullopt;
    const Eigen::Vector3d intersection_point = ray.start_pos + ray_scaling * ray.dir;

    const bool inside = normal.dot( ( tri.v2 - tri.v1 ).cross( intersection_point - tri.v1 ) ) >= 0 and
                        normal.dot( ( tri.v3 - tri.v2 ).cross( intersection_point - tri.v2 ) ) >= 0 and
                        normal.dot( ( tri.v1 - tri.v3 ).cross( intersection_point - tri.v3 ) ) >= 0;

    return inside ? std::optional<Eigen::Vector3d>( intersection_point ) : std::nullopt;
}

std::optional<TracePoint> traceRayOnTet( const topology::TetMeshCombinatorialMap& map,
                                         const topology::Cell& start_cell,
                                         const Ray<3>& ray,
                                         const std::vector<Normal>& normals )
{
    const auto face_ids = indexingOrError( map, 2 );
    // The face on the input dart is the location that we start from.
    // Check all the other faces for intersection with the ray.
    LOG( LOG_TRACING ) << "Tracing on tet " << start_cell.dart().id() << " from ray " << ray.start_pos.transpose()
                       << " -> " << ray.dir.transpose() << std::endl;

    // Avoid attempting intersection with faces adjacent to start_cell.
    topology::LocalCellMarker m( 2 );
    iterateAdjacentCells( map, start_cell, 2, [&]( const auto& I ) { m.mark( map, I ); return true; } );
    bool found_it = false;
    std::pair<topology::Face, Eigen::Vector3d> out;
    iterateAdjacentCells( map, topology::Volume( start_cell.dart() ), 2, [&]( const topology::Face& adj_face ) {
        if( m.isMarked( adj_face ) ) return true;
        LOG( LOG_TRACING ) << "| Face " << adj_face.dart().id() << std::endl;
        const Triangle<3> tri = triangleOfFace( map, adj_face );
        LOG( LOG_TRACING ) << "| | Triangle<3>: " << tri.v1.transpose() << " | " << tri.v2.transpose() << " | "
                           << tri.v3.transpose() << std::endl;

        const std::optional<Eigen::Vector3d> intersection =
            intersectionOf( ray, tri, normals.at( face_ids( adj_face ) ).get( adj_face.dart() ) );

        if( intersection.has_value() )
        {
            LOG( LOG_TRACING ) << "| | Found intersection: " << intersection.value().transpose() << std::endl;
            const auto maybe_phi3 = phi( map, 3, adj_face.dart() );
            out.first = topology::Face( maybe_phi3.has_value() ? maybe_phi3.value() : adj_face.dart() );
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

std::optional<TracePoint> traceCellAverageField( const topology::TetMeshCombinatorialMap& map,
                                                 const topology::Cell& c,
                                                 const Eigen::Vector3d& start_point,
                                                 const Eigen::Matrix3Xd& field,
                                                 const std::vector<Normal>& normals )
{
    const auto volume_ids = indexingOrError( map, 3 );
    const auto face_ids = indexingOrError( map, 2 );
    LOG( LOG_TRACING ) << "| | | Tracing on cell average field\n";
    Eigen::Vector3d ave_field = Eigen::Vector3d::Zero();
    iterateAdjacentCells( map, c, 3, [&]( const topology::Volume& v ) {
        ave_field += field.col( volume_ids( v ) );
        return true;
    } );

    return tracingStartCell(
               map,
               c,
               [&]( const topology::Face& f ) { return normals.at( face_ids( f ) ).get( f.dart() ); },
               [&]( const topology::Volume& ) { return ave_field; } )
        .and_then( [&]( const topology::Cell& start_cell ) {
            return traceRayOnTet( map, start_cell, Ray<3>{ start_point, ave_field }, normals );
        } );
}

SimplicialComplex traceField( const topology::TetMeshCombinatorialMap& map,
                              const topology::Cell& start_cell,
                              const Eigen::Vector3d& start_point,
                              const Eigen::Matrix3Xd& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output )
{
    const auto volume_ids = indexingOrError( map, 3 );
    const auto face_ids = indexingOrError( map, 2 );

    topology::Cell curr_cell;
    Eigen::Vector3d curr_point = start_point;

    SimplicialComplex complex;
    SimplicialComplex debug_tets;
    complex.points.push_back( curr_point );

    // Figure out the start face. If there is a face adjacent to this cell that works, use that.
    if( start_cell.dim() == 2 )
    {
        curr_cell = start_cell;
    }
    else
    {
        const auto normals_func = [&]( const topology::Face& f ) { return normals.at( face_ids( f ) ).get( f.dart() ); };
        const auto grads_func = [&]( const topology::Volume& v ) { return field.col( volume_ids( v ) ); };
        const auto edge_func = [&]( const topology::Face& f ) {
            const Triangle<3> tri = triangleOfFace( map, f );
            return ( tri.v3 + tri.v2 ) * 0.5 - tri.v1;
        };

        const std::optional<topology::Cell> adjusted_start_cell =
            tracingStartCell( map, start_cell, normals_func, grads_func ).or_else( [&]() {
                return tracingAverageStartCell( map, start_cell, normals_func, grads_func, edge_func );
            } );

        if( adjusted_start_cell.has_value() )
        {
            curr_cell = adjusted_start_cell.value();
        }
        else throw( std::runtime_error( "Untraceable field at start" ) );
    }

    size_t n = 0;
    LOG( LOG_TRACING ) << "Starting a trace\n";
    do
    {
        const topology::Volume curr_vol( curr_cell.dart() );
        const Ray<3> search_ray( { curr_point, field.col( volume_ids( curr_vol ) ) } );
        if( debug_output ) interiorTracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n );
        std::optional<TracePoint> next_point = traceRayOnTet( map, curr_cell, search_ray, normals );
        if( not next_point.has_value() )
        {
            next_point = traceCellAverageField( map, curr_cell, curr_point, field, normals );
            if( not next_point.has_value() )
            {
                throw( std::runtime_error( "Untraceable field at step " + std::to_string( n ) ) );
            }
        }
        curr_cell = next_point.value().first;
        curr_point = next_point.value().second;
        complex.points.push_back( curr_point );
        complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
        n++;
        if( n > 2000 ) throw std::runtime_error( "Never ending tracing loop" );
    } while( not onBoundary( map, curr_cell.dart() ) );

    return complex;
}


////////////////////////////////////////
////////////////////////////////////////
/// BOUNDARY TRACING
////////////////////////////////////////

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
    return barycentricIntersectionOf( ray, line ).transform( [&]( const double& u ) -> Eigen::Vector2d {
        return ( 1 - u ) * line.start_pos + u * line.end_pos;
    } );
}

/// Tri holds the vertex positions at start_cell.dart(), its phi1, and its phi(1,1), in that order.
std::optional<std::pair<topology::Edge, double>> run2dIntersectionOnTriangle( const topology::CombinatorialMap& map,
                                                                              const Triangle<2>& tri,
                                                                              const topology::Cell& start_cell,
                                                                              const Ray<2>& ray )
{
    const auto segment = [&]( const size_t i ) -> Segment<2> {
        switch( i )
        {
            case 0: return Segment<2>( { tri.v1, tri.v2 } );
            case 1: return Segment<2>( { tri.v2, tri.v3 } );
            case 2: return Segment<2>( { tri.v3, tri.v1 } );
            default: throw std::runtime_error( "Why did you put in the wrong number?" );
        }
    };
    const auto run_intersection = [&]( const topology::Dart& d_next, const Segment<2>& line ) {
        return barycentricIntersectionOf( ray, line )
            .transform( [&]( const double& u ) -> std::pair<topology::Edge, double> {
                const auto maybe_next_face_dart = phi( map, 2, d_next );
                const topology::Edge e( maybe_next_face_dart.value_or( d_next ) );
                return { e, maybe_next_face_dart.has_value() ? 1.0 - u : u };
            } );
    };

    // Avoid attempting intersection with edges adjacent to start_cell.
    topology::LocalCellMarker m( 1 );
    iterateAdjacentCells( map, start_cell, 1, [&]( const auto& d ) { m.mark( map, d ); return true; } );

    topology::Dart iter_dart = start_cell.dart();
    for( size_t i = 0; i < 3; i++, iter_dart = phi( map, 1, iter_dart ).value() )
    {
        if( not m.isMarked( topology::Cell( iter_dart, 1 ) ) )
        {
            const auto possible_out = run_intersection( iter_dart, segment( i ) );
            if( possible_out.has_value() ) return possible_out;
        }
    }
    return std::nullopt;
}

Eigen::Vector2d gradient2d( const Triangle<2> tri, const Eigen::Vector3d& field, const double twice_area )
{
    const auto segment = [&]( const size_t i ) -> Segment<2> {
        switch( i )
        {
            case 0: return Segment<2>( { tri.v1, tri.v2 } );
            case 1: return Segment<2>( { tri.v2, tri.v3 } );
            case 2: return Segment<2>( { tri.v3, tri.v1 } );
            default: throw std::runtime_error( "Why did you put in the wrong number?" );
        }
    };
    const Eigen::Rotation2Dd rot90( std::numbers::pi / 2.0 );

    const auto grad_s_i = [&]( const Segment<2>& edge_i ) -> Eigen::Vector2d {
        const auto edge_diff = edge_i.end_pos - edge_i.start_pos;
        return 1.0 / twice_area * ( rot90 * edge_diff );
    };

    const Eigen::Vector2d gradient =
        field( 0 ) * grad_s_i( segment( 1 ) ) + field( 1 ) * grad_s_i( segment( 2 ) ) + field( 2 ) * grad_s_i( segment( 0 ) );

    return gradient;
}

std::optional<std::pair<topology::Edge, double>>
    traceGradientOnTri( const topology::CombinatorialMap& map,
                        const VertexPositionsFunc& positions,
                        const topology::Cell& start_cell,
                        const double edge_barycentric_coord,
                        const Eigen::VectorXd& field_values )
{
    const topology::Face f( start_cell.dart() );
    const Triangle<3> tri3d = triangleOfFace<3>( map, positions, f );

    const topology::Dart& d = start_cell.dart();

    // move the triangle into the xy plane
    const Eigen::Vector3d e0 = ( tri3d.v2 - tri3d.v1 ).normalized();
    const Eigen::Vector3d e1 = ( tri3d.v3 - tri3d.v1 - e0.dot( tri3d.v3 - tri3d.v1 ) * e0 ).normalized();

    const Eigen::Vector2d v0_2d( 0, 0 );
    const Eigen::Vector2d v1_2d( ( tri3d.v2 - tri3d.v1 ).norm(), 0 );
    const Eigen::Vector2d v2_2d( e0.dot( tri3d.v3 - tri3d.v1 ), e1.dot( tri3d.v3 - tri3d.v1 ) );

    const double twice_area = v1_2d( 0 ) * v2_2d( 1 );

    const auto vertex_ids = indexingOrError( map, 0 );

    // calculate the gradient in the xy plane
    const Eigen::Vector3d face_field = field_values( { vertex_ids( topology::Vertex( d ) ),
                                                       vertex_ids( topology::Vertex( phi( map, 1, d ).value() ) ),
                                                       vertex_ids( topology::Vertex( phi( map, -1, d ).value() ) ) } );

    const Eigen::Vector2d gradient = gradient2d( { v0_2d, v1_2d, v2_2d }, face_field, twice_area );

    // iterate the edges and perform line ray intersections
    const Eigen::Vector2d start_coord = start_cell.dim() == 0 ? Eigen::Vector2d::Zero() : Eigen::Vector2d( edge_barycentric_coord * v1_2d );
    const Ray<2> ray( { start_coord, gradient } );

    return run2dIntersectionOnTriangle( map, {v0_2d, v1_2d, v2_2d }, start_cell, ray );
}

std::optional<std::pair<topology::Edge, double>>
    traceGradientOnTriPair( const topology::CombinatorialMap& map,
                            const VertexPositionsFunc& positions,
                            const topology::Cell& start_cell,
                            const double edge_barycentric_coord,
                            const Eigen::VectorXd& field_values )
{
    return phi( map, { 2, -1 }, start_cell.dart() ).and_then( [&]( const topology::Dart& opp_d ) {
        const topology::Face f( start_cell.dart() );
        const Triangle<3> tri3d = triangleOfFace<3>( map, positions, f );
        const Eigen::Vector3d opp_v = positions( topology::Vertex( opp_d ) );

        // calculate the gradient in the xy plane
        const topology::Dart& d = start_cell.dart();

        // move the triangle into the xy plane
        const Eigen::Vector3d e0 = ( tri3d.v2 - tri3d.v1 ).normalized();
        const Eigen::Vector3d e1 = ( tri3d.v3 - tri3d.v1 - e0.dot( tri3d.v3 - tri3d.v1 ) * e0 ).normalized();
        const Eigen::Vector3d e1_prime = -( opp_v - tri3d.v1 - e0.dot( opp_v - tri3d.v1 ) * e0 ).normalized();

        const Eigen::Vector2d v0_2d( 0, 0 );
        const Eigen::Vector2d v1_2d( ( tri3d.v2 - tri3d.v1 ).norm(), 0 );
        const Eigen::Vector2d v2_2d( e0.dot( tri3d.v3 - tri3d.v1 ), e1.dot( tri3d.v3 - tri3d.v1 ) );
        const Eigen::Vector2d v3_2d( e0.dot( opp_v - tri3d.v1 ), e1_prime.dot( opp_v - tri3d.v1 ) );

        const auto vertex_ids = indexingOrError( map, 0 );
        // calculate the gradient in the xy plane
        const Eigen::Vector4d faces_field =
            field_values( { vertex_ids( topology::Vertex( d ) ),
                            vertex_ids( topology::Vertex( phi( map, 1, d ).value() ) ),
                            vertex_ids( topology::Vertex( phi( map, -1, d ).value() ) ),
                            vertex_ids( topology::Vertex( opp_d ) ) } );

        const double twice_area_0 = v1_2d( 0 ) * v2_2d( 1 );
        const double twice_area_1 = -v1_2d( 0 ) * v3_2d( 1 );

        const Eigen::Vector2d gradient_0 = gradient2d( { v0_2d, v1_2d, v2_2d }, faces_field( { 0, 1, 2 } ), twice_area_0 );
        const Eigen::Vector2d gradient_1 = gradient2d( { v0_2d, v3_2d, v1_2d }, faces_field( { 0, 3, 1 } ), twice_area_1 );

        const Eigen::Vector2d gradient = ( gradient_0 + gradient_1 ) * 0.5;

        // iterate the edges and perform line ray intersections
        const Eigen::Vector2d start_coord =
            start_cell.dim() == 0 ? Eigen::Vector2d::Zero() : Eigen::Vector2d( edge_barycentric_coord * v1_2d );
        const Ray<2> ray( { start_coord, gradient } );

        if( gradient( 1 ) >= 0 )
        {
            return run2dIntersectionOnTriangle( map, {v0_2d, v1_2d, v2_2d }, start_cell, ray );
        }
        else
        {
            if( start_cell.dim() == 1 )
            {
                const topology::Cell start_cell_opp( phi( map, 1, opp_d ).value(), start_cell.dim() );
                return run2dIntersectionOnTriangle( map, {v1_2d, v0_2d, v3_2d }, start_cell_opp, ray );
            }
            else if( start_cell.dim() == 0 )
            {
                const topology::Cell start_cell_opp( phi( map, -1, opp_d ).value(), start_cell.dim() );
                return run2dIntersectionOnTriangle( map, {v0_2d, v3_2d, v1_2d }, start_cell_opp, ray );
            }
            else throw std::runtime_error( "Can't trace from a face on a face" );
        }
    } );
}

Trace traceBoundaryField( const topology::CombinatorialMap& map,
                          const topology::Cell& start_cell,
                          const double& start_point,
                          const Eigen::VectorXd& field,
                          const VertexPositionsFunc& positions,
                          const bool debug_output )
{
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto expand_barycentric = [&map, &positions]( const topology::Cell& c, const double coord ) -> Eigen::Vector3d {
        if( c.dim() == 0 ) return positions( c );
        else if( c.dim() == 1 )
        {
            const topology::Dart& d = c.dart();
            const Eigen::Vector3d pos1 = positions( topology::Vertex( d ) );
            const Eigen::Vector3d pos2 = positions( topology::Vertex( phi( map, 1, d ).value() ) );

            return ( 1.0 - coord ) * pos1 + coord * pos2;
        }
        else throw std::runtime_error( "Don't pass triangle or higher as the start cell" );
    };
    const auto expand_field = [&map, &field, &vertex_ids]( const topology::Cell& c, const double coord ) -> double {
        if( c.dim() == 0 ) return field( vertex_ids( c ) );
        else if( c.dim() == 1 )
        {
            const topology::Dart& d = c.dart();
            const double pos1 = field( vertex_ids( topology::Vertex( d ) ) );
            const double pos2 = field( vertex_ids( topology::Vertex( phi( map, 1, d ).value() ) ) );

            return ( 1.0 - coord ) * pos1 + coord * pos2;
        }
        else throw std::runtime_error( "Don't pass triangle or higher as the start cell" );
    };

    topology::Cell curr_cell;
    double curr_point = start_point;

    SimplicialComplex traced_line;
    SimplicialComplex debug_tris;
    std::vector<double> debug_field_values; // This is at the tris, not the trace
    std::vector<double> field_values;       // At the verts of the trace
    std::vector<topology::Cell> tracing_faces;

    if( start_cell.dim() == 1 )
    {
        curr_cell = start_cell;
    }
    else
    {
        const auto normals_func = [&]( const topology::Edge& e ) -> Eigen::Vector3d {
            // Rotate edge vector 90 degrees about the normal axis
            const Eigen::Vector3d edge_vector =
                positions( topology::Vertex( phi( map, 1, e.dart() ).value() ) ) - positions( topology::Vertex( e.dart() ) );

            const Eigen::Vector3d tri_normal = triangleNormal( triangleOfFace<3>( map, positions, topology::Face( e.dart() ) ) );

            return Eigen::AngleAxis<double>( std::numbers::pi / 2.0, tri_normal ) * edge_vector;
        };
        const auto edge_func = [&]( const topology::Edge& e ) -> Eigen::Vector3d {
            return positions( topology::Vertex( phi( map, 1, e.dart() ).value() ) ) - positions( topology::Vertex( e.dart() ) );
        };
        const auto grads_func = [&]( const topology::Face& f ) -> Eigen::Vector3d {
            const topology::Dart& d = f.dart();
            const auto face_field = field( { vertex_ids( topology::Vertex( d ) ),
                                             vertex_ids( topology::Vertex( phi( map, 1, d ).value() ) ),
                                             vertex_ids( topology::Vertex( phi( map, -1, d ).value() ) ) } );
            return gradient( triangleOfFace<3>( map, positions, f ), face_field );
        };
        const std::optional<topology::Cell> adjusted_start_cell =
            tracingStartCell( map, start_cell, normals_func, grads_func ).or_else( [&]() {
                return tracingAverageStartCell( map, start_cell, normals_func, grads_func, edge_func );
            } );

        if( adjusted_start_cell.has_value() )
        {
            curr_cell = adjusted_start_cell.value();
        }
        else throw( std::runtime_error( "Untraceable field at start" ) );
    }

    size_t n = 0;
    traced_line.points.push_back( expand_barycentric( start_cell, start_point ) );
    field_values.push_back( expand_field( start_cell, start_point ) );
    do
    {
        LOG( LOG_TRACING ) << n << " - Attempting trace from " << curr_cell << ", " << curr_point << std::endl;
        const topology::Face curr_face( curr_cell.dart() );
        if( debug_output )
            boundaryTracingDebugOutput( map, positions, curr_cell, curr_point, field, traced_line, debug_tris, debug_field_values, n );
        std::optional<std::pair<topology::Edge, double>> next_point =
            traceGradientOnTri( map, positions, curr_cell, curr_point, field );

        if( not next_point.has_value() )
        {
            LOG( LOG_TRACING ) << "Didn't get a tracing match.\n";
            next_point = traceGradientOnTriPair( map, positions, curr_cell, curr_point, field );
            if( not next_point.has_value() )
            {
                throw std::runtime_error( "No intersection on boundary tracing n=" + std::to_string( n ) );
                break;
            }
        }

        curr_cell = next_point.value().first;
        curr_point = next_point.value().second;
        traced_line.points.push_back( expand_barycentric( curr_cell, curr_point ) );
        field_values.push_back( expand_field( curr_cell, curr_point ) );
        traced_line.simplices.push_back( { traced_line.points.size() - 2, traced_line.points.size() - 1 } );
        tracing_faces.push_back( curr_face );
        n++;
    } while( not onBoundary( map, curr_cell.dart() ) and n < 2000 );
    // FIXME: This number being the maximum trace length is completely arbitrary, and I assume at some point this will
    // break an actual tracing that is supposed to be this long.
    if( n > 1900 )
        throw std::runtime_error( "Unending tracing loop Cell(" + std::to_string( start_cell.dart().id() ) + ", " +
                                  std::to_string( start_cell.dim() ) + ")" );

    return { traced_line, field_values, tracing_faces };
}

}