#pragma once
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include<set>
#include<map>
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

    static std::optional<double> edgeCrossingPoint( const cgogn::CMap3& map,
                                                    const cgogn::CMap3::Edge& e,
                                                    const double desired_value,
                                                    const std::function<double( const VertexId& )>& vertex_vals )
    {
        const double val1 = vertex_vals( cgogn::index_of( map, cgogn::CMap3::Vertex( e.dart_ ) ) );
        const double val2 = vertex_vals( cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, e.dart_ ) ) ) );

        if( not ( val1 < desired_value and val2 > desired_value ) and not ( val1 > desired_value and val2 < desired_value ) )
            return {};
 
        return ( desired_value - val1 ) / ( val2 - val1 );
    }

    static std::pair< std::vector<Simplex>, std::vector<Eigen::Vector3d> > isosurface(
        const cgogn::CMap3& map,
        const double desired_value,
        const std::function<double( const VertexId& )>& vertex_vals )
    {
        std::vector<Simplex> simplices_out;
        std::vector<Eigen::Vector3d> points_out;

        //auto edge_crossings = cgogn::add_attribute<VertexId, cgogn::CMap3::Vertex>( map, "__new_vertices" );
        std::map<uint32_t, VertexId> edge_crossings;

        const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
        cgogn::foreach_cell( map, [&]( cgogn::CMap3::Edge e ) {
            const std::optional<double> crossing_point = edgeCrossingPoint( map, e, desired_value, vertex_vals );

            if( crossing_point )
            {
                edge_crossings.emplace( cgogn::index_of( map, e ), VertexId( points_out.size() ) );

                const auto& d = e.dart_;
                const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
                const Eigen::Vector3d& pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );

                const auto& t = crossing_point.value();
                points_out.push_back( ( 1 - t ) * pos1 + t * pos2 );
            }
            return true;
        } );

        /*
          1. Check every edge to see if it crosses the iso surface, and mark those that do. (LOOK INTO CGOGN PROPERTIES (not worth it))
          These edges also need a new vertex created at the crossing, stored as a vertex id in the new complex and a 3d position. (DONE ABOVE)

          2. Check every tet to see if it has three adjacent edges that cross the surface. (std::map::findOptional would be useful here)
          When this happens add a new triangle to the new complex using the points stored on those edges.
        */

       return { simplices_out, points_out };
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