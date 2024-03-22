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

constexpr bool LOG_TRACING = 0;

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
/// INTERIOR TRACING
////////////////////////////////////////

std::optional<Eigen::Vector3d> intersectionOf( const Ray<3>& ray,
                                               const Triangle<3>& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal )
{
    const Eigen::Vector3d normal = maybe_normal.value_or( triangleNormal( tri ) );
    const double ray_scaling = normal.dot( tri.v1 - ray.start_pos ) / normal.dot( ray.dir );
    LOG( LOG_TRACING ) << "| | Ray scaling: " << ray_scaling << std::endl;
    // TODO: Is there a better approach here? If it is within that tolerance it should be an edge I guess?
    if( ray_scaling < 1e-13 ) return std::nullopt;
    const Eigen::Vector3d intersection_point = ray.start_pos + ray_scaling * ray.dir;

    const bool inside = normal.dot( ( tri.v2 - tri.v1 ).cross( intersection_point - tri.v1 ) ) >= 0 and
                        normal.dot( ( tri.v3 - tri.v2 ).cross( intersection_point - tri.v2 ) ) >= 0 and
                        normal.dot( ( tri.v1 - tri.v3 ).cross( intersection_point - tri.v3 ) ) >= 0;

    return inside ? std::optional<Eigen::Vector3d>( intersection_point ) : std::nullopt;
}

std::optional<TracePoint> traceRayOnTet( const topology::TetMeshCombinatorialMap& map,
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

static constexpr bool LOG_TRACING_CELL = false;
std::optional<topology::Cell> tracingStartCell( const topology::CombinatorialMap& map,
                                                const topology::Cell& lower_dim_cell,
                                                const std::function<Eigen::Vector3d( const topology::Cell& )>& normal,
                                                const std::function<Eigen::Vector3d( const topology::Cell& )>& grad )
{
    LOG( LOG_TRACING_CELL ) << "Searching for tracing cell...\n";
    const size_t dim = map.dim();

    LOG( LOG_TRACING_CELL ) << "Cell dim: " << lower_dim_cell.dim() << std::endl;
    std::optional<topology::Cell> out;
    iterateAdjacentCells( map, lower_dim_cell, dim, [&]( const topology::Cell& elem ) {
        const Eigen::Vector3d this_grad = grad( elem );
        LOG( LOG_TRACING_CELL ) << "Grad: " << this_grad.transpose() << std::endl;
        const bool found_cell =
            iterateAdjacentCellsOfRestrictedCell( map,
                                                  topology::Cell( elem.dart(), lower_dim_cell.dim() ),
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
            out.emplace( topology::Cell( elem.dart(), dim - 1 ) );
            return false;
        }
        return true;
    } );
    return out;
}

std::optional<TracePoint> traceCellAverageField( const topology::TetMeshCombinatorialMap& map,
                                                 const topology::Cell& c,
                                                 const Eigen::Vector3d& start_point,
                                                 const Eigen::MatrixX3d& field,
                                                 const std::vector<Normal>& normals )
{
    LOG( LOG_TRACING ) << "| | | Tracing on cell average field\n";
    Eigen::Vector3d ave_field = Eigen::Vector3d::Zero();
    iterateAdjacentCells( map, c, 3, [&]( const topology::Volume& v ) {
        ave_field += field.row( map.elementId( v ) ).transpose();
        return true;
    } );

    return tracingStartCell(
               map,
               c,
               [&]( const topology::Face& f ) { return normals.at( map.faceId( f ) ).get( f.dart() ); },
               [&]( const topology::Volume& ) { return ave_field; } )
        .and_then( [&]( const topology::Face& start_f ) {
            return traceRayOnTet( map, topology::Volume( start_f.dart() ), Ray<3>{ start_point, ave_field }, normals );
        } );
}

SimplicialComplex traceField( const topology::TetMeshCombinatorialMap& map,
                              const topology::Cell& start_cell,
                              const Eigen::Vector3d& start_point,
                              const Eigen::MatrixX3d& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output )
{
    topology::Face curr_face;
    Eigen::Vector3d curr_point = start_point;

    SimplicialComplex complex;
    SimplicialComplex debug_tets;
    complex.points.push_back( curr_point );

    // Figure out the start face. If there is a face adjacent to this cell that works, use that.
    if( start_cell.dim() == 2 )
    {
        curr_face = start_cell;
    }
    else
    {
        const std::optional<topology::Face> start_face = tracingStartCell(
            map,
            start_cell,
            [&]( const topology::Face& f ) { return normals.at( map.faceId( f ) ).get( f.dart() ); },
            [&]( const topology::Volume& v ) { return field.row( map.elementId( v ) ).transpose(); } );

        if( start_face.has_value() )
        {
            curr_face = start_face.value();
        }
        else
        {
            // If there was not a face adjacent to start_cell that works, trace the cell-average field once.
            std::optional<TracePoint> second_point = traceCellAverageField( map, start_cell, curr_point, field, normals );
            if( not second_point.has_value() )
            {
                throw( std::runtime_error( "Untraceable field at start" ) );
            }
            curr_face = second_point.value().first;
            curr_point = second_point.value().second;
            complex.points.push_back( curr_point );
            complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
        }
    }

    size_t n = 0;
    LOG( LOG_TRACING ) << "Starting a trace\n";
    do
    {
        const topology::Volume curr_vol( curr_face.dart() );
        const Ray<3> search_ray( { curr_point, field.row( map.elementId( curr_vol ) ) } );
        if( debug_output ) tracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n++ );
        std::optional<TracePoint> next_point = traceRayOnTet( map, curr_vol, search_ray, normals );
        if( not next_point.has_value() )
        {
            next_point = traceCellAverageField( map, topology::Face( curr_vol.dart() ), curr_point, field, normals );
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


////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
/// BOUNDARY TRACING

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


struct TriPairTraceResult
{
    bool on_forward_face;
    bool on_forward_edge;
    double edge_barycoord;
};

std::optional<TriPairTraceResult> traceGradientOnTriPair( const Triangle<3>& tri3d,
                                                          const Eigen::Vector3d& opp_v,
                                                          const double edge_barycentric_coord,
                                                          const Eigen::Ref<const Eigen::Vector4d> field )
{
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

    const auto run_intersection = [&]( const bool f, const bool i, const Segment<2>& line ) {
        return barycentricIntersectionOf( ray, line )
            .and_then( [&]( const double& u ) -> std::optional<TriPairTraceResult> {
                return TriPairTraceResult{ f, i, u };
            } );
    };

    if( gradient( 1 ) >= 0 )
    {
        return run_intersection( FORWARD, FORWARD, { v2_2d, v3_2d } )
            .or_else( [&]() {
                return run_intersection( FORWARD, BACKWARD, { v3_2d, v1_2d } );
            } );
    }
    else
    {
        return run_intersection( BACKWARD, FORWARD, { v4_2d, v2_2d } )
            .or_else( [&]() {
                return run_intersection( BACKWARD, BACKWARD, { v1_2d, v4_2d } );
            } );
    }
}

void tracingDebugOutput( const topology::CombinatorialMap& map,
                         const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& positions,
                         const topology::Face& curr_face,
                         const double start_pos,
                         const Eigen::Ref<const Eigen::VectorXd> field,
                         const SimplicialComplex& line,
                         SimplicialComplex& tris,
                         std::vector<double>& vertex_values,
                         const size_t n )
{
    const Triangle<3> tri3d = triangleOfFace( map, positions, curr_face );

    const auto& d = curr_face.dart();
    const auto field_values = field( { map.vertexId( topology::Vertex( d ) ).id(),
                                       map.vertexId( topology::Vertex( phi( map, 1, d ).value() ) ).id(),
                                       map.vertexId( topology::Vertex( phi( map, -1, d ).value() ) ).id() } );

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
    const auto grad = gradient( tri3d, field_values );
    const Ray<3> ray( { ( 1.0 - start_pos ) * tri3d.v1 + start_pos * tri3d.v2, grad } );
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

std::optional<std::pair<topology::Edge, double>>
    traceGradientOnTri( const topology::CombinatorialMap& map,
                        const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& positions,
                        const topology::Face& f,
                        const double edge_barycentric_coord,
                        const Eigen::VectorXd& field_values )
{
    const Triangle<3> tri3d = triangleOfFace( map, positions, f );

    const topology::Dart& d = f.dart();
    const auto field = field_values( { map.vertexId( topology::Vertex( d ) ).id(),
                                       map.vertexId( topology::Vertex( phi( map, 1, d ).value() ) ).id(),
                                       map.vertexId( topology::Vertex( phi( map, -1, d ).value() ) ).id() } );

    return traceGradientOnTri( tri3d, edge_barycentric_coord, field )
        .and_then( [&]( const std::pair<bool, double>& pr ) {
            const int first_op = pr.first ? 1 : -1;
            const auto maybe_next_face_dart = phi( map, { first_op, 2 }, f.dart() );
            const topology::Edge e( maybe_next_face_dart.value_or( phi( map, first_op, f.dart() ).value() ) );
            return std::optional<std::pair<topology::Edge, double>>( { e, maybe_next_face_dart.has_value() ? 1.0 - pr.second : pr.second } );
        } );
}

std::optional<std::pair<topology::Edge, double>>
    traceGradientOnTriPair( const topology::CombinatorialMap& map,
                            const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& positions,
                            const topology::Face& f,
                            const double edge_barycentric_coord,
                            const Eigen::VectorXd& field_values )
{
    return phi( map, { 2, -1 }, f.dart() ).and_then( [&]( const topology::Dart& opp_d ) {
        const Triangle<3> tri3d = triangleOfFace( map, positions, f );
        const Eigen::Vector3d& opp_v = positions( topology::Vertex( opp_d ) );

        // calculate the gradient in the xy plane
        const topology::Dart& d = f.dart();
        const auto field = field_values( { map.vertexId( topology::Vertex( d ) ).id(),
                                        map.vertexId( topology::Vertex( phi( map, 1, d ).value() ) ).id(),
                                        map.vertexId( topology::Vertex( phi( map, -1, d ).value() ) ).id(),
                                        map.vertexId( topology::Vertex( phi( map, { 2, -1 }, d ).value() ) ).id() } );

        return traceGradientOnTriPair( tri3d, opp_v, edge_barycentric_coord, field )
            .and_then( [&]( const TriPairTraceResult& result ) {
                if( result.on_forward_face )
                {
                    const int first_op = result.on_forward_edge ? 1 : -1;
                    const auto maybe_next_face_dart = phi( map, { first_op, 2 }, f.dart() );
                    const topology::Edge e( maybe_next_face_dart.value_or( phi( map, first_op, f.dart() ).value() ) );
                    return std::optional<std::pair<topology::Edge, double>>( { e, 1.0 - result.edge_barycoord } );
                }
                else
                {
                    const int second_op = result.on_forward_edge ? 1 : -1;
                    const auto maybe_next_face_dart = phi( map, { 2, second_op, 2 }, f.dart() );
                    const topology::Edge e( maybe_next_face_dart.value_or( phi( map, { 2, second_op }, f.dart() ).value() ) );
                    return std::optional<std::pair<topology::Edge, double>>( { e, 1.0 - result.edge_barycoord } );
                }
            } );
    } );
}

SimplicialComplex traceBoundaryField( const topology::CombinatorialMap& map,
                                      const topology::Edge& e,
                                      const double& start_point,
                                      const Eigen::VectorXd& field,
                                      const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& positions,
                                      const bool debug_output,
                                      const std::function<void( const topology::Face& )>& face_callback )
{
    const auto expand_barycentric = [&map, &positions]( const topology::Edge& e, const double coord ) {
        const topology::Dart& d = e.dart();
        const Eigen::Vector3d& pos1 = positions( topology::Vertex( d ) );
        const Eigen::Vector3d& pos2 = positions( topology::Vertex( phi( map, 1, d ).value() ) );

        return ( 1.0 - coord ) * pos1 + coord * pos2;
    };

    topology::Edge curr_edge = e;
    double curr_point = start_point;

    SimplicialComplex traced_line;
    SimplicialComplex debug_tris;
    std::vector<double> debug_field_values;
    size_t n = 0;

    traced_line.points.push_back( expand_barycentric( curr_edge, curr_point ) );
    do
    {
        const topology::Face curr_face( curr_edge.dart() );
        face_callback( curr_face );
        if( debug_output )
            tracingDebugOutput( map, positions, curr_face, curr_point, field, traced_line, debug_tris, debug_field_values, n++ );
        std::optional<std::pair<topology::Edge, double>> next_point =
            traceGradientOnTri( map, positions, curr_face, curr_point, field );

        if( not next_point.has_value() )
        {
            next_point = traceGradientOnTriPair( map, positions, curr_face, curr_point, field );
            if( not next_point.has_value() )
            {
                throw std::runtime_error( "No intersection on boundary tracing" );
                break;
            }
        }

        curr_edge = next_point.value().first;
        curr_point = next_point.value().second;
        traced_line.points.push_back( expand_barycentric( curr_edge, curr_point ) );
        traced_line.simplices.push_back( { traced_line.points.size() - 2, traced_line.points.size() - 1 } );
    } while( not onBoundary( map, curr_edge.dart() ) );

    return traced_line;
}