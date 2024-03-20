#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <set>
#include <iomanip>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <SweepInput.hpp>
#include <Logging.hpp>

Eigen::Vector3d triangleNormal( const Triangle<3>& tri )
{
    return ( tri.v2 - tri.v1 ).cross( tri.v3 - tri.v1 ).normalized();
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

Eigen::Vector3d triangleNormal( const topology::TetMeshCombinatorialMap& map, const topology::Face& f )
{
    return triangleNormal( triangleOfFace( map, f ) );
}

Eigen::Vector3d centroid( const Triangle<3>& tri )
{
    return ( tri.v1 + tri.v2 + tri.v3 ) / 3;
}

Eigen::Vector3d centroid( const topology::TetMeshCombinatorialMap& map, const topology::Face& f )
{
    return centroid( triangleOfFace( map, f ) );
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

double dihedralCotangent( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals )
{
    const topology::Dart& d1 = e.dart();
    const topology::Dart d2 = phi( map, 2, e.dart() ).value();
    const Eigen::Vector3d n1 = normals.at( map.faceId( topology::Face( d1 ) ) ).get( d1 );
    const Eigen::Vector3d n2 = normals.at( map.faceId( topology::Face( d2 ) ) ).reversed( d2 );

    const double cos_theta = n1.dot( n2 );
    return cos_theta / std::sqrt( 1 - cos_theta * cos_theta );
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

Eigen::Vector3d expandBarycentric( const topology::CombinatorialMap& map,
                                   const VertexPositionsFunc& positions,
                                   const topology::Cell& start_face,
                                   const BarycentricPoint& coord )
{
    Eigen::Vector3d out = Eigen::Vector3d::Zero();
    const auto& simplex_verts = coord.simplex.vertices();
    iterateAdjacentCells( map, start_face, 0, [&]( const topology::Vertex& v ) {
        const VertexId vid = map.vertexId( v );
        const size_t idx = std::distance( simplex_verts.begin(), std::find( simplex_verts.begin(), simplex_verts.end(), vid ) );
        if( idx >= simplex_verts.size() )
        {
            std::cout << coord.simplex << std::endl;
            std::cout << vid << std::endl;
            std::cout << idx << std::endl;
            throw std::runtime_error( "Face and simplex don't match" );
        }
        out += coord.point( idx ) * positions( v );
        return true;
    } );
    return out;
}