#pragma once
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include<set>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <SweepInput.hpp>
#include<Logging.hpp>

class SimplexUtilities
{
    public:
    static Eigen::Vector3d triangleNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3 )
    {
        return ( v2 - v1 ).cross( v3 - v1 ).normalized();
    }

    static Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
    {
        const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
        const cgogn::Dart& d = f.dart_;
        const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
        const Eigen::Vector3d& pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
        const Eigen::Vector3d& pos3 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
        return triangleNormal( pos1, pos2, pos3 );
    }

    static double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
    {
        const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
        const cgogn::Dart& d = e.dart_;
        const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
        const Eigen::Vector3d& pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
        return ( pos2 - pos1 ).norm();
    }

    static double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
    {
        const Eigen::Vector3d n1 = triangleNormal( map, cgogn::CMap3::Face( e.dart_ ) );
        const Eigen::Vector3d n2 = triangleNormal( map, cgogn::CMap3::Face( cgogn::phi2( map, e.dart_ ) ) );

        const double cos_theta = abs( n1.dot( n2 ) );
        return cos_theta / std::sqrt( 1 - cos_theta * cos_theta );
    }

    // FIXME: This doesn't belong here
    static void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map )
    {
        cgogn::io::VolumeImportData import;
        import.reserve( sweep_input.points.size(), sweep_input.simplices.size() );

        for( const auto& tet : sweep_input.simplices )
        {
            import.volumes_types_.push_back( cgogn::io::VolumeType::Tetra );
            for( size_t i = 0; i < 4; i++ )
            {
                import.volumes_vertex_indices_.push_back( tet.vertex( i ).id() );
            }
        }

        import.vertex_position_ = sweep_input.points;

        import_volume_data( map, import );

        cgogn::index_cells<cgogn::CMap3::Edge>( map );
    }
};