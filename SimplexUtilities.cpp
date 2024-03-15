#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <set>
#include <iomanip>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
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

Triangle<3> triangleOfFace( const topology::TetMeshCombinatorialMap& map, const topology::Face& f )
{
    const SimplicialComplex& complex = map.simplicialComplex();
    const auto vertex_position = [&]( const topology::Vertex& v ) -> const Eigen::Vector3d& {
        return complex.points.at( map.vertexId( v ).id() );
    };

    return triangleOfFace( map, vertex_position, f );
}

Triangle<3> triangleOfFace( const topology::CombinatorialMap& map,
                            const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& vertex_position,
                            const topology::Face& f )
{
    const topology::Dart& d = f.dart();

    const Eigen::Vector3d& pos1 = vertex_position( topology::Vertex( d ) );
    const Eigen::Vector3d& pos2 = vertex_position( topology::Vertex( phi( map, 1, d ).value() ) );
    const Eigen::Vector3d& pos3 = vertex_position( topology::Vertex( phi( map, -1, d ).value() ) );

    return Triangle<3>{ pos1, pos2, pos3 };
}

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    return triangleNormal( triangleOfFace( map, f ) );
}

Eigen::Vector3d triangleNormal( const topology::TetMeshCombinatorialMap& map, const topology::Face& f )
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

Eigen::Vector3d centroid( const topology::TetMeshCombinatorialMap& map, const topology::Face& f )
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

std::vector<Normal> faceNormals( const topology::TetMeshCombinatorialMap& map )
{
    const size_t n_faces = cellCount( map, 2 );
    std::vector<Normal> normals( n_faces );

    iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
        const auto fid = map.faceId( f );
        normals[fid] = Normal( map, f.dart(), triangleNormal( map, f ) );
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

double edgeLength( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e )
{
    const SimplicialComplex& complex = map.simplicialComplex();
    const topology::Dart& d = e.dart();
    const auto vertex_position = [&]( const topology::Vertex& v ) {
        return complex.points.at( map.vertexId( v ).id() );
    };
    const Eigen::Vector3d& pos1 = vertex_position( topology::Vertex( d ) );
    const Eigen::Vector3d& pos2 = vertex_position( topology::Vertex( phi( map, 1, d ).value() ) );
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

double dihedralCotangent( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals )
{
    const topology::Dart& d1 = e.dart();
    const topology::Dart d2 = phi( map, 2, e.dart() ).value();
    const Eigen::Vector3d n1 = normals.at( map.faceId( topology::Face( d1 ) ) ).get( d1 );
    const Eigen::Vector3d n2 = normals.at( map.faceId( topology::Face( d2 ) ) ).reversed( d2 );

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

Eigen::Vector3d gradient( const topology::TetMeshCombinatorialMap& map,
                          const topology::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals )
{
    using namespace topology;
    const SimplicialComplex& complex = map.simplicialComplex();
    const auto vertex_position = [&]( const topology::Vertex& v ) {
        return complex.points.at( map.vertexId( v ).id() );
    };
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    iterateAdjacentCells( map, v, 2, [&]( const Face& f ) {
        const Vertex op_vert( phi( map, {2, -1}, f.dart() ).value() );
        const VertexId op_vert_id = map.vertexId( op_vert );
        const Eigen::Vector3d& op_vert_pos = vertex_position( op_vert );
        const Eigen::Vector3d& face_vert_pos = vertex_position( Vertex( f.dart() ) );
        const auto fid = map.faceId( f );
        const Eigen::Vector3d& normal = normals.at( fid ).get( f.dart() );
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

Eigen::MatrixX3d gradients( const topology::TetMeshCombinatorialMap& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals )
{
    Eigen::MatrixX3d result( cellCount( map, 3 ), 3 );
    iterateCellsWhile( map, 3, [&]( const topology::Volume& v ) {
        result.row( map.elementId( v ) ) = gradient( map, v, field_values, normals ).transpose();
        return true;
    } );

    return result;
}

Eigen::Vector3d gradient( const Triangle<3>& tri3d, const Eigen::Ref<const Eigen::Vector3d> field_values )
{
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

    const Eigen::Vector3d grad_3d = gradient( 0 ) * e1 + gradient( 1 ) * e2;
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