#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <set>
#include <iomanip>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <SweepInput.hpp>
#include <Logging.hpp>

Eigen::Vector3d triangleNormal( const Triangle<3>& tri )
{
    return ( tri.v2 - tri.v1 ).cross( tri.v3 - tri.v1 ).normalized();
}

Triangle<3> triangleOfFace( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = f.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    const Eigen::Vector3d& pos3 =
        cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );

    return Triangle<3>{ pos1, pos2, pos3 };
}

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    return triangleNormal( triangleOfFace( map, f ) );
}

Eigen::Vector3d centroid( const Triangle<3>& tri )
{
    return ( tri.v1 + tri.v2 + tri.v3 ) / 3;
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

Eigen::Vector3d gradient( const Triangle<3>& tri3d,
                          const Eigen::Ref<const Eigen::Vector3d> field_values )
{
    // std::cout << "Triangle: " << tri3d.v1.transpose() << " | " << tri3d.v2.transpose() << " | " << tri3d.v3.transpose() << std::endl;
    // std::cout << field_values.transpose() << std::endl;
    // Basis for the plane of the triangle
    const Eigen::Vector3d e1 = ( tri3d.v2 - tri3d.v1 ).normalized();
    const Eigen::Vector3d e2 = ( tri3d.v3 - tri3d.v1 - e1.dot( tri3d.v3 - tri3d.v1 ) * e1 ).normalized();
    // std::cout << std::setprecision(10) << "e1: " << e1.transpose() << std::endl << "e2: " << e2.transpose() << std::endl;

    const Eigen::Vector2d v1_2d( 0, 0 );
    const Eigen::Vector2d v2_2d( ( tri3d.v2 - tri3d.v1 ).norm(), 0 );
    const Eigen::Vector2d v3_2d( e1.dot( tri3d.v3 - tri3d.v1 ), e2.dot( tri3d.v3 - tri3d.v1 ) );

    // calculate the gradient in the xy plane
    const double& f1 = field_values( 0 );
    const double& f2 = field_values( 1 );
    const double& f3 = field_values( 2 );

    const Eigen::Rotation2Dd rot90( std::numbers::pi / 2.0 );

    // std::cout << Eigen::Matrix2d( rot90 ) << std::endl << std::endl;
    const double twice_area = v2_2d( 0 ) * v3_2d( 1 );

    const auto grad_s_i = [&]( const Segment<2>& edge_i ) -> Eigen::Vector2d {
        const auto edge_diff = edge_i.end_pos - edge_i.start_pos;
        return 1.0 / twice_area * ( rot90 * edge_diff );
    };

    // std::cout << "Edge 1: " << grad_s_i( { v1_2d, v2_2d } ).transpose() << std::endl;
    // std::cout << "Edge 2: " << grad_s_i( { v2_2d, v3_2d } ).transpose() << std::endl;
    // std::cout << "Edge 3: " << grad_s_i( { v3_2d, v1_2d } ).transpose() << std::endl;

    const Eigen::Vector2d gradient =
        f1 * grad_s_i( { v2_2d, v3_2d } ) + f2 * grad_s_i( { v3_2d, v1_2d } ) + f3 * grad_s_i( { v1_2d, v2_2d } );

    // std::cout << "2d Gradient: " << gradient.transpose() << std::endl;

    const Eigen::Vector3d grad_3d = gradient( 0 ) * e1 + gradient( 1 ) * e2;
    // std::cout << "Gradient: " << grad_3d.transpose() << std::endl << std::endl;
    return grad_3d;
}

void addTriangleNoDuplicateChecking( SimplicialComplex& complex, const Triangle<3>& tri )
{
    const size_t offset = complex.points.size();
    complex.points.push_back( tri.v1 );
    complex.points.push_back( tri.v2 );
    complex.points.push_back( tri.v3 );
    complex.simplices.push_back( Simplex( offset, offset + 1, offset + 2 ) );
}

void mapFromInput( const SimplicialComplex& mesh, cgogn::CMap3& map )
{
    cgogn::io::VolumeImportData import;
    import.reserve( mesh.points.size(), mesh.simplices.size() );

    for( const auto& tet : mesh.simplices )
    {
        import.volumes_types_.push_back( cgogn::io::VolumeType::Tetra );
        for( size_t i = 0; i < 4; i++ )
        {
            import.volumes_vertex_indices_.push_back( tet.vertex( i ).id() );
        }
    }

    import.vertex_position_ = mesh.points;

    import_volume_data( map, import );

    cgogn::index_cells<cgogn::CMap3::Edge>( map );
    cgogn::index_cells<cgogn::CMap3::Face>( map );
    cgogn::index_cells<cgogn::CMap3::Volume>( map );
}