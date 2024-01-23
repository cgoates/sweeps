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

struct Triangle
{
    const Eigen::Ref<const Eigen::Vector3d> v1;
    const Eigen::Ref<const Eigen::Vector3d> v2;
    const Eigen::Ref<const Eigen::Vector3d> v3;
};

Triangle triangleOfFace( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

Eigen::Vector3d triangleNormal( const Triangle& tri );

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

Eigen::Vector3d centroid( const Triangle& tri );

Eigen::Vector3d centroid( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

std::vector<Normal> faceNormals( const cgogn::CMap3& map );

double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals );

Eigen::Vector3d gradient( const cgogn::CMap3& map,
                          const cgogn::CMap3::Volume& v,
                          const Eigen::VectorXd& field_values,
                          const std::vector<Normal>& normals );

// FIXME: This doesn't belong here
void mapFromInput( const SimplicialComplex& mesh, cgogn::CMap3& map );