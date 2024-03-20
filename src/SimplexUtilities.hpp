#pragma once
#include <CombinatorialMapMethods.hpp>
#include <SweepInput.hpp>

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

template<unsigned int DIM>
struct Triangle
{
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> v1;
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> v2;
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> v3;
};

template <unsigned int DIM> struct Segment
{
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> start_pos;
    const Eigen::Ref<const Eigen::Matrix<double, DIM, 1>> end_pos;
};

void addTriangleNoDuplicateChecking( SimplicialComplex& complex, const Triangle<3>& tri );

Triangle<3> triangleOfFace( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );
Triangle<3> triangleOfFace( const topology::CombinatorialMap& map,
                            const std::function<const Eigen::Vector3d&( const topology::Vertex& )>& vertex_position,
                            const topology::Face& f );

Eigen::Vector3d triangleNormal( const Triangle<3>& tri );

Eigen::Vector3d triangleNormal( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );

Eigen::Vector3d centroid( const Triangle<3>& tri );

Eigen::Vector3d centroid( const topology::TetMeshCombinatorialMap& map, const topology::Face& f );

std::vector<Normal> faceNormals( const topology::TetMeshCombinatorialMap& map );

double edgeLength( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e );

double dihedralCotangent( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const topology::TetMeshCombinatorialMap& map,
                          const topology::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals );

Eigen::MatrixX3d gradients( const topology::TetMeshCombinatorialMap& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const Triangle<3>& tri3d, const Eigen::Ref<const Eigen::Vector3d> field_values );

using VertexPositionsFunc = std::function<const Eigen::Vector3d&( const topology::Vertex& )>;
Eigen::Vector3d expandBarycentric( const topology::CombinatorialMap& map,
                                   const VertexPositionsFunc& positions,
                                   const topology::Cell& start_face,
                                   const BarycentricPoint& coord );