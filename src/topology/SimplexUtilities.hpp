#pragma once
#include <CombinatorialMapMethods.hpp>
#include <SweepInput.hpp>
#include <VertexPositionsFunc.hpp>

namespace topology
{
    class TetMeshCombinatorialMap;
}

class Normal
{
    public:
    Normal( const topology::CombinatorialMap& map, const topology::Dart& dart, const Eigen::Vector3d& normal )
        : mAlignedDarts( { dart, phi( map, 1, dart ).value(), phi( map, -1, dart ).value() } ),
          mNormal( normal )
    {}
    Normal() : mAlignedDarts( { 0, 0, 0 } ) {}

    Eigen::Vector3d get( const topology::Dart& dart ) const
    {
        if( std::find( mAlignedDarts.begin(), mAlignedDarts.end(), dart ) == mAlignedDarts.end() )
        {
            return -1 * mNormal;
        }
        return mNormal;
    }
    Eigen::Vector3d reversed( const topology::Dart& dart ) const
    {
        if( std::find( mAlignedDarts.begin(), mAlignedDarts.end(), dart ) == mAlignedDarts.end() )
        {
            return mNormal;
        }
        return -1 * mNormal;
    }

    private:
    std::array<topology::Dart, 3> mAlignedDarts;
    Eigen::Vector3d mNormal;
};

struct Tetrahedron
{
    const Eigen::Vector3d v1;
    const Eigen::Vector3d v2;
    const Eigen::Vector3d v3;
    const Eigen::Vector3d v4;
};

template<int DIM>
struct Triangle
{
    const Eigen::Matrix<double, DIM, 1> v1;
    const Eigen::Matrix<double, DIM, 1> v2;
    const Eigen::Matrix<double, DIM, 1> v3;
};

template <int DIM> struct Segment
{
    const Eigen::Matrix<double, DIM, 1> start_pos;
    const Eigen::Matrix<double, DIM, 1> end_pos;
};

void addEdgeNoDuplicateChecking( SimplicialComplex& complex,
                                 const topology::CombinatorialMap& map,
                                 const VertexPositionsFunc& pos,
                                 const topology::Edge& e );

void addTriangleNoDuplicateChecking( SimplicialComplex& complex, const Triangle<3>& tri );
void addAllTriangles( SimplicialComplex& complex, const topology::CombinatorialMap& cmap, const VertexPositionsFunc& pos );

void addTetNoDuplicateChecking( SimplicialComplex& complex,
                                const topology::CombinatorialMap& map,
                                const VertexPositionsFunc& pos,
                                const topology::Volume& vol );

Triangle<3> triangleOfFace( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );

template <int DIM>
Triangle<DIM> triangleOfFace( const topology::CombinatorialMap& map,
                              const VertexPositionsFunc& vertex_position,
                              const topology::Face& f );

Tetrahedron tetOfVolume( const topology::CombinatorialMap& map, const VertexPositionsFunc& v_positions, const topology::Volume& v );

Eigen::Vector3d triangleNormal( const Triangle<3>& tri );

Eigen::Vector3d triangleNormal( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );

Eigen::Vector3d centroid( const Triangle<3>& tri );
Eigen::Vector3d centroid( const Tetrahedron& tet );

Eigen::Vector3d centroid( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );

std::vector<Normal> faceNormals( const topology::TetMeshCombinatorialMap& map );

double edgeLength( const topology::CombinatorialMap& map, const VertexPositionsFunc& positions, const topology::Edge& e );
double tetVolume( const Tetrahedron& tet );

double dihedralCosine( const topology::CombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals );
double dihedralCotangent( const topology::CombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals );

Eigen::Vector3d circumcenter( const Tetrahedron& tet );
template<int DIM> Eigen::Matrix<double, DIM, 1> circumcenter( const Triangle<DIM>& tri );

Eigen::Vector3d gradient( const topology::TetMeshCombinatorialMap& map,
                          const topology::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals );

Eigen::Matrix3Xd gradients( const topology::TetMeshCombinatorialMap& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals );

Eigen::Matrix3Xd gradientsWithBoundaryCorrection( const topology::TetMeshCombinatorialMap& map,
                                                  const topology::CombinatorialMap& sides,
                                                  const Eigen::VectorXd& field_values,
                                                  const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const Triangle<3>& tri3d, const Eigen::Ref<const Eigen::Vector3d> field_values );

Eigen::Vector3d expandBarycentric( const topology::CombinatorialMap& map,
                                   const VertexPositionsFunc& positions,
                                   const topology::Cell& start_face,
                                   const BarycentricPoint& coord );

bool isInverted( const topology::CombinatorialMap& map,
                 const topology::Volume& v,
                 const VertexPositionsFunc& positions );

std::optional<Eigen::Vector3d> invertTriangleMap( const Triangle<2>& tri, const Eigen::Vector2d& point );
std::optional<Eigen::Vector3d> invertTriangleMap( const Triangle<3>& tri, const Eigen::Vector3d& point );

/// @brief Returns the closest barycentric coordinates in a triangle to the input point, as well as the distance from
/// the triangle to the point if the point is not in the triangle.
/// @param tri   The triangle to invert
/// @param point The point
/// @return The barycentric coordinates of the closest point in the triangle to point, as well as the distance if nonzero.
std::pair<Eigen::Vector3d, std::optional<double>> invertTriangleMapOrClosestPoint( const Triangle<3>& tri, const Eigen::Vector3d& point );
std::pair<Eigen::Vector3d, std::optional<double>> invertTriangleMapOrClosestPoint( const Triangle<2>& tri, const Eigen::Vector2d& point );

// Find { 1 - t, t } for pt = ( 1 - t ) * pt0 + t * pt1
std::pair<double, double> inverseLinear( const double pt0, const double pt1, const double pt );