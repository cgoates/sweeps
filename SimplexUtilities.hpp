#pragma once
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/types/cell_marker.h>
#include <SweepInput.hpp>

struct SweepInput;

class Normal
{
    public:
    Normal( const cgogn::CMap3& map, const cgogn::Dart& dart, const Eigen::Vector3d& normal ) :
        mOppositeDarts( { cgogn::phi3( map, dart ), cgogn::phi<3, 1>( map, dart ), cgogn::phi<3, -1>( map, dart ) } ),
        mNormal( normal )
    {}
    Normal() {}

    Eigen::Vector3d get( const cgogn::Dart& dart ) const
    {
        if( std::find( mOppositeDarts.begin(), mOppositeDarts.end(), dart ) != mOppositeDarts.end() ) return -1 * mNormal;
        return mNormal;
    }

    private:
    std::array<cgogn::Dart, 3> mOppositeDarts;
    Eigen::Vector3d mNormal;
};
class SimplexUtilities
{
    public:
    static Eigen::Vector3d triangleNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3 );

    static Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

    static std::vector<Normal> faceNormals( const cgogn::CMap3& map );

    static double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

    static double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e, const std::vector<Normal>& normals );

    static Eigen::Vector3d gradient( const cgogn::CMap3& map,
                                     const cgogn::CMap3::Volume& v,
                                     const std::function<double(const cgogn::CMap3::Vertex&)>& field_values,
                                     const std::function<const Eigen::Vector3d&(const cgogn::CMap3::Face&)>& inward_normals );

    // FIXME: This doesn't belong here
    static void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map );
};