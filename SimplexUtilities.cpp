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

Eigen::Vector3d triangleNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3 )
{
    return ( v2 - v1 ).cross( v3 - v1 ).normalized();
}

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = f.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    const Eigen::Vector3d& pos3 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
    return triangleNormal( pos1, pos2, pos3 );
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
        gradient += field_values( op_vert_id.id() ) * normal.dot( op_vert_pos - face_vert_pos ) * normal;
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

void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map )
{
    cgogn::io::VolumeImportData import;
    import.reserve( sweep_input.points.size(), sweep_input.simplices.size() );

    for( const auto& tet : sweep_input.simplices )
    {
        import.volumes_types_.push_back( cgogn::io::VolumeType::Tetra );
        for( size_t i = 0; i < 4; i++ )
        {
            import.volumes_vertex_indices_.push_back( tet.vertex( i ).id() );
        }
    }

    import.vertex_position_ = sweep_input.points;

    import_volume_data( map, import );

    cgogn::index_cells<cgogn::CMap3::Edge>( map );
    cgogn::index_cells<cgogn::CMap3::Face>( map );
}