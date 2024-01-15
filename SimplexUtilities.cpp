#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include <set>
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

constexpr bool LOG_TRACING = 1;

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
}