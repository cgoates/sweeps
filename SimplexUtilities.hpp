#pragma once
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/types/cell_marker.h>
#include <SweepInput.hpp>

class Normal
{
    public:
    Normal( const cgogn::CMap3& map, const cgogn::Dart& dart, const Eigen::Vector3d& normal )
        : mOppositeDarts( { cgogn::phi3( map, dart ), cgogn::phi<3, 1>( map, dart ), cgogn::phi<3, -1>( map, dart ) } ),
          mNormal( normal )
    {}
    Normal() {}

    Eigen::Vector3d get( const cgogn::Dart& dart ) const
    {
        if( std::find( mOppositeDarts.begin(), mOppositeDarts.end(), dart ) != mOppositeDarts.end() )
            return -1 * mNormal;
        return mNormal;
    }

    private:
    std::array<cgogn::Dart, 3> mOppositeDarts;
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

Triangle<3> triangleOfFace( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

Eigen::Vector3d triangleNormal( const Triangle<3>& tri );

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

Eigen::Vector3d centroid( const Triangle<3>& tri );

Eigen::Vector3d centroid( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

std::vector<Normal> faceNormals( const cgogn::CMap3& map );

double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const cgogn::CMap3& map,
                          const cgogn::CMap3::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals );

Eigen::MatrixX3d gradients( const cgogn::CMap3& map,
                            const Eigen::VectorXd& field_values,
                            const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const Triangle<3>& tri3d, const Eigen::Ref<const Eigen::Vector3d> field_values );

// FIXME: This doesn't belong here
void mapFromInput( const SimplicialComplex& mesh, cgogn::CMap3& map );