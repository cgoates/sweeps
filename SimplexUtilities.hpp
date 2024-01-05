#pragma once
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <SweepInput.hpp>

struct SweepInput;

class SimplexUtilities
{
    public:
    static Eigen::Vector3d triangleNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3 );

    static Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f );

    static double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

    static double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

    static Eigen::Vector3d gradient( const cgogn::CMap3& map,
                                     const cgogn::CMap3::Volume& v,
                                     const std::function<double(const cgogn::CMap3::Vertex&)>& field_values,
                                     const std::function<const Eigen::Vector3d&(const cgogn::CMap3::Face&)>& inward_normals );

    // FIXME: This doesn't belong here
    static void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map );
};