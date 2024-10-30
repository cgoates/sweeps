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
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto vertex_position = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
        return complex.points.at( vertex_ids( v ) );
    };

    return triangleOfFace<3>( map, vertex_position, f );
}

template<unsigned int DIM>
Triangle<DIM> triangleOfFace( const topology::CombinatorialMap& map,
                              const VertexPositionsFunc& vertex_position,
                              const topology::Face& f )
{
    const topology::Dart& d = f.dart();

    const Eigen::Matrix<double, DIM, 1> pos1 = vertex_position( topology::Vertex( d ) ).head<DIM>();
    const Eigen::Matrix<double, DIM, 1> pos2 = vertex_position( topology::Vertex( phi( map, 1, d ).value() ) ).head<DIM>();
    const Eigen::Matrix<double, DIM, 1> pos3 = vertex_position( topology::Vertex( phi( map, -1, d ).value() ) ).head<DIM>();

    return Triangle<DIM>{ pos1, pos2, pos3 };
}
template
Triangle<3> triangleOfFace( const topology::CombinatorialMap&,
                            const VertexPositionsFunc&,
                            const topology::Face& );
template
Triangle<2> triangleOfFace( const topology::CombinatorialMap&,
                            const VertexPositionsFunc&,
                            const topology::Face& );

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
    const auto face_ids = indexingOrError( map, 2 );
    const size_t n_faces = cellCount( map, 2 );
    std::vector<Normal> normals( n_faces );

    iterateCellsWhile( map, 2, [&]( const topology::Face& f ) {
        const auto fid = face_ids( f );
        normals[fid] = Normal( map, f.dart(), triangleNormal( map, f ) );
        return true;
    } );
    return normals;
}

double edgeLength( const topology::CombinatorialMap& map, const VertexPositionsFunc& vertex_position, const topology::Edge& e )
{
    const topology::Dart& d = e.dart();
    const Eigen::Vector3d pos1 = vertex_position( topology::Vertex( d ) );
    const Eigen::Vector3d pos2 = vertex_position( topology::Vertex( phi( map, 1, d ).value() ) );
    return ( pos2 - pos1 ).norm();
}

double dihedralCotangent( const topology::CombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals )
{
    const auto face_ids = indexingOrError( map, 2 );
    const topology::Dart& d1 = e.dart();
    const topology::Dart d2 = phi( map, 2, e.dart() ).value();
    const Eigen::Vector3d n1 = normals.at( face_ids( topology::Face( d1 ) ) ).get( d1 );
    const Eigen::Vector3d n2 = normals.at( face_ids( topology::Face( d2 ) ) ).reversed( d2 );

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
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto face_ids = indexingOrError( map, 2 );
    const auto vertex_position = [&]( const topology::Vertex& v ) {
        return complex.points.at( vertex_ids( v ) );
    };
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    iterateAdjacentCells( map, v, 2, [&]( const Face& f ) {
        const Vertex op_vert( phi( map, {2, -1}, f.dart() ).value() );
        const VertexId op_vert_id = vertex_ids( op_vert );
        const Eigen::Vector3d op_vert_pos = vertex_position( op_vert );
        const Eigen::Vector3d face_vert_pos = vertex_position( Vertex( f.dart() ) );
        const auto fid = face_ids( f );
        const Eigen::Vector3d& normal = normals.at( fid ).get( f.dart() );
        gradient += field_values( op_vert_id.id() ) * normal / normal.dot( op_vert_pos - face_vert_pos );
        return true;
    } );

    return gradient;
}

Eigen::Matrix3Xd gradients( const topology::TetMeshCombinatorialMap& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals )
{
    const auto volume_ids = indexingOrError( map, 3 );
    Eigen::Matrix3Xd result( 3, cellCount( map, 3 ) );
    iterateCellsWhile( map, 3, [&]( const topology::Volume& v ) {
        result.col( volume_ids( v ) ) = gradient( map, v, field_values, normals );
        return true;
    } );

    return result;
}

Eigen::Matrix3Xd gradientsWithBoundaryCorrection( const topology::TetMeshCombinatorialMap& map,
                                                  const topology::CombinatorialMap& sides,
                                                  const Eigen::VectorXd& field_values,
                                                  const std::vector<Normal>& normals )
{
    const auto volume_ids = indexingOrError( map, 3 );
    const auto face_ids = indexingOrError( map, 2 );
    Eigen::Matrix3Xd result = gradients( map, field_values, normals );
    iterateCellsWhile( sides, 2, [&]( const topology::Face& f ) {
        const size_t col = volume_ids( topology::Volume( f.dart() ) );
        const size_t fid = face_ids( f );
        const Eigen::Vector3d outward_normal = normals.at( fid ).reversed( f.dart() );
        const auto normal_portion = result.col( col ).dot( outward_normal );
        result.col( col ) -= normal_portion * outward_normal;
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

void addEdgeNoDuplicateChecking( SimplicialComplex& complex,
                                 const topology::CombinatorialMap& map,
                                 const VertexPositionsFunc& pos,
                                 const topology::Edge& e )
{
    const size_t offset = complex.points.size();
    complex.points.push_back( pos( e.dart() ) );
    complex.points.push_back( pos( phi( map, 1, e.dart() ).value() ) );
    complex.simplices.push_back( Simplex( offset, offset + 1 ) );
}

void addTriangleNoDuplicateChecking( SimplicialComplex& complex, const Triangle<3>& tri )
{
    const size_t offset = complex.points.size();
    complex.points.push_back( tri.v1 );
    complex.points.push_back( tri.v2 );
    complex.points.push_back( tri.v3 );
    complex.simplices.push_back( Simplex( offset, offset + 1, offset + 2 ) );
}

void addTetNoDuplicateChecking( SimplicialComplex& complex,
                                const topology::CombinatorialMap& map,
                                const VertexPositionsFunc& pos,
                                const topology::Volume& vol )
{
    const size_t offset = complex.points.size();
    complex.points.push_back( pos( vol.dart() ) );
    complex.points.push_back( pos( phi( map, 1, vol.dart() ).value() ) );
    complex.points.push_back( pos( phi( map, -1, vol.dart() ).value() ) );
    complex.points.push_back( pos( phi( map, {2, -1}, vol.dart() ).value() ) );
    complex.simplices.push_back( Simplex( offset, offset + 1, offset + 2, offset + 3 ) );
}

Eigen::Vector3d expandBarycentric( const topology::CombinatorialMap& map,
                                   const VertexPositionsFunc& positions,
                                   const topology::Cell& start_face,
                                   const BarycentricPoint& coord )
{
    Eigen::Vector3d out = Eigen::Vector3d::Zero();
    const auto vertex_ids = indexingOrError( map, 0 );
    const auto& simplex_verts = coord.simplex.vertices();
    iterateAdjacentCells( map, start_face, 0, [&]( const topology::Vertex& v ) {
        const VertexId vid = vertex_ids( v );
        const size_t idx = std::distance( simplex_verts.begin(), std::find( simplex_verts.begin(), simplex_verts.end(), vid ) );
        if( idx >= simplex_verts.size() )
        {
            std::cerr << coord.simplex << std::endl;
            std::cerr << vid << std::endl;
            std::cerr << idx << std::endl;
            throw std::runtime_error( "Face and simplex don't match" );
        }
        out += coord.point( idx ) * positions( v );
        return true;
    } );
    return out;
}

bool isInverted( const topology::CombinatorialMap& map,
                 const topology::Volume& v,
                 const VertexPositionsFunc& positions )
{
    const auto tri = triangleOfFace<3>( map, positions, topology::Face( v.dart() ) );
    const auto normal = triangleNormal( tri );
    const auto inversion_question = ( positions( topology::Vertex( phi( map, { 2, -1 }, v.dart() ).value() ) ) - tri.v1 ).dot( normal );
    if( inversion_question < 0 )
    {
        std::cout << "How inverted? " << inversion_question << std::endl;
        return true;
    }
    return false;
}

std::optional<Eigen::Vector3d> invertTriangleMap( const Triangle<2>& tri, const Eigen::Vector2d& point )
{
    // See https://gamedev.stackexchange.com/a/63203
    const Eigen::Vector2d diff0 = tri.v2 - tri.v1;
    const Eigen::Vector2d diff1 = tri.v3 - tri.v1;
    const Eigen::Vector2d diff_point = point - tri.v1;
    const auto perpendicular_dot_product = []( const Eigen::Vector2d& a, const Eigen::Vector2d& b ){
        return a( 0 ) * b( 1 ) - a( 1 ) * b( 0 );
    };
    const double denominator = 1.0 / perpendicular_dot_product( diff0, diff1 );
    Eigen::Vector3d out = Eigen::Vector3d::Zero();
    out( 1 ) = denominator * perpendicular_dot_product( diff_point, diff1 );
    if( out( 1 ) > 1.0 or out( 1 ) < -1e-15 ) return std::nullopt;
    out( 2 ) = denominator * perpendicular_dot_product( diff0, diff_point );
    if( out( 2 ) > 1.0 or out( 2 ) < -1e-15 ) return std::nullopt;
    out( 0 ) = 1 - ( out( 1 ) + out( 2 ) );
    if( out( 0 ) > 1.0 or out( 0 ) < -1e-15 ) return std::nullopt;

    return out;
}

std::pair<double, double> inverseLinear( const double pt0, const double pt1, const double pt )
{
    const double t = ( pt - pt0 ) / ( pt1 - pt0 );
    return { 1 - t, t };
}