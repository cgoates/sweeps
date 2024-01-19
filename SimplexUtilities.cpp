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
#include <SweepInput.hpp>
#include <Logging.hpp>

Eigen::Vector3d triangleNormal( const Triangle& tri )
{
    return ( tri.v2 - tri.v1 ).cross( tri.v3 - tri.v1 ).normalized();
}

Triangle triangleOfFace( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = f.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    const Eigen::Vector3d& pos3 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );

    return Triangle{ pos1, pos2, pos3 };
}

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    return triangleNormal( triangleOfFace( map, f ) );
}

Eigen::Vector3d centroid( const Triangle& tri )
{
    return 0.33 * ( tri.v1 + tri.v2 + tri.v3 );
}

Eigen::Vector3d centroid( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    return centroid( triangleOfFace( map, f ) );
}

std::vector<Normal> faceNormals( const cgogn::CMap3& map )
{
    const size_t n_faces = cgogn::nb_cells<cgogn::CMap3::Face>( map );
    std::vector<Normal> normals( n_faces );

    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
        const auto fid = cgogn::index_of( map, f );
        normals[fid] = Normal( map, f.dart_, triangleNormal( map, f ) );
        return true;
    } );
    return normals;
}

double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = e.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    return ( pos2 - pos1 ).norm();
}

double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals )
{
    const auto get_normal = [&]( cgogn::CMap3::Face f ) {
        const auto fid = cgogn::index_of( map, f );
        return normals.at( fid ).get( f.dart_ );
    };
    const Eigen::Vector3d n1 = get_normal( cgogn::CMap3::Face( e.dart_ ) );
    const Eigen::Vector3d n2 = get_normal( cgogn::CMap3::Face( cgogn::phi<2, 3>( map, e.dart_ ) ) );

    const double cos_theta = n1.dot( n2 );
    return cos_theta / std::sqrt( 1 - cos_theta * cos_theta );
}

Eigen::Vector3d gradient( const cgogn::CMap3& map,
                          const cgogn::CMap3::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals )
{
    using Vertex = cgogn::CMap3::Vertex;
    const auto position = cgogn::get_attribute<Eigen::Vector3d, Vertex>( map, "position" );
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    cgogn::foreach_incident_face( map, v, [&]( cgogn::CMap3::Face f ) {
        const Vertex op_vert( cgogn::phi<2, -1>( map, f.dart_ ) );
        const VertexId op_vert_id( cgogn::index_of( map, op_vert ) );
        const Eigen::Vector3d& op_vert_pos = cgogn::value<Eigen::Vector3d>( map, position, op_vert );
        const Eigen::Vector3d& face_vert_pos = cgogn::value<Eigen::Vector3d>( map, position, Vertex( f.dart_ ) );
        const auto fid = cgogn::index_of( map, f );
        const Eigen::Vector3d& normal = normals.at( fid ).get( f.dart_ );
        gradient += field_values( op_vert_id.id() ) * normal / normal.dot( op_vert_pos - face_vert_pos );
        return true;
    } );

    return gradient;
}

Eigen::MatrixX3d gradients( const cgogn::CMap3& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals )
{
    Eigen::MatrixX3d result( cgogn::nb_cells<cgogn::CMap3::Volume>( map ), 3 );
    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Volume v ) {
        result.row( cgogn::index_of( map, v ) ) = gradient( map, v, field_values, normals ).transpose();
        return true;
    } );

    return result;
}

constexpr bool LOG_TRACING = 0;

std::optional<Eigen::Vector3d> intersectionOf( const Ray& ray,
                                               const Triangle& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal )
{
    const Eigen::Vector3d normal = maybe_normal.value_or( triangleNormal( tri ) );
    const double ray_scaling = normal.dot( tri.v1 - ray.start_pos ) / normal.dot( ray.dir );
    LOG( LOG_TRACING ) << "| | Ray scaling: " << ray_scaling << std::endl;
    if( ray_scaling <= 0 ) return std::nullopt;
    const Eigen::Vector3d intersection_point =
        ray.start_pos + ray_scaling * ray.dir;

    const bool inside = normal.dot( ( tri.v2 - tri.v1 ).cross( intersection_point - tri.v1 ) ) >= 0 and
                        normal.dot( ( tri.v3 - tri.v2 ).cross( intersection_point - tri.v2 ) ) >= 0 and
                        normal.dot( ( tri.v1 - tri.v3 ).cross( intersection_point - tri.v3 ) ) >= 0;

    return inside ? std::optional<Eigen::Vector3d>( intersection_point ) : std::nullopt;
}

std::optional<TracePoint> traceRayOnTet( const cgogn::CMap3& map,
                                         const cgogn::CMap3::Volume& v,
                                         const Ray& ray,
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
        const Triangle tri = triangleOfFace( map, adj_face );
        LOG( LOG_TRACING ) << "| | Triangle: " << tri.v1.transpose() << " | " << tri.v2.transpose() << " | "
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
    if( not found_it ) return std::nullopt;
    return out;
}

void tracingDebugOutput( const cgogn::CMap3& map,
                         const cgogn::CMap3::Volume& v,
                         const Ray& ray,
                         const SimplicialComplex& line,
                         SimplicialComplex& tets,
                         const size_t n )
{
    // Add the triangles of the current tet to tets
    foreach_incident_face( map, v, [&]( cgogn::CMap3::Face f ) {
        const Triangle from_face = triangleOfFace( map, f );
        const size_t next_vid = tets.points.size();
        tets.simplices.push_back( { next_vid + 0, next_vid + 1, next_vid + 2 } );
        tets.points.push_back( from_face.v1 );
        tets.points.push_back( from_face.v2 );
        tets.points.push_back( from_face.v3 );
        return true;
    } );

    // Create a simplicial complex with just one face: the triangle it's coming from
    const Triangle from_face = triangleOfFace( map, cgogn::CMap3::Face( v.dart_ ) );
    const SimplicialComplex from_face_complex( { { {0, 1, 2} }, { from_face.v1, from_face.v2, from_face.v3 } } );
    const io::VTKOutputObject from_face_output( from_face_complex );

    // Create a simplicial complex for the ray
    const SimplicialComplex ray_complex( { { {0} }, { ray.start_pos } } );
    io::VTKOutputObject ray_output( ray_complex );
    ray_output.addVertexField( "ray", ray.dir.transpose() );

    const io::VTKOutputObject line_output( line );
    const io::VTKOutputObject tets_output( tets );

    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << n;
    std::string n_str(ss.str());

    io::outputSimplicialFieldToVTK( from_face_output, "from_face_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( ray_output, "ray_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( line_output, "line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( tets_output, "tets_" + n_str + ".vtu" );
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
        const Ray search_ray( {curr_point, field.row( index_of( map, curr_vol ) )} );
        if( debug_output ) tracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n++ );
        const std::optional<TracePoint> next_point = traceRayOnTet( map, curr_vol, search_ray, normals );
        if( not next_point.has_value() ) break;
        curr_face = next_point.value().first;
        curr_point = next_point.value().second;
        complex.points.push_back( curr_point );
        complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
    }

    return complex;
}

void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map )
{
    cgogn::io::VolumeImportData import;
    import.reserve( sweep_input.mesh.points.size(), sweep_input.mesh.simplices.size() );

    for( const auto& tet : sweep_input.mesh.simplices )
    {
        import.volumes_types_.push_back( cgogn::io::VolumeType::Tetra );
        for( size_t i = 0; i < 4; i++ )
        {
            import.volumes_vertex_indices_.push_back( tet.vertex( i ).id() );
        }
    }

    import.vertex_position_ = sweep_input.mesh.points;

    import_volume_data( map, import );

    cgogn::index_cells<cgogn::CMap3::Edge>( map );
    cgogn::index_cells<cgogn::CMap3::Face>( map );
    cgogn::index_cells<cgogn::CMap3::Volume>( map );
}