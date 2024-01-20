#include <Tracing.hpp>
#include <SimplexUtilities.hpp>
#include <iostream>
#include <Simplex.hpp>
#include <VTKOutput.hpp>
#include <set>
#include <iomanip>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <SweepInput.hpp>
#include <Logging.hpp>

constexpr bool LOG_TRACING = 0;

std::optional<Eigen::Vector3d> intersectionOf( const Ray& ray,
                                               const Triangle& tri,
                                               std::optional<const Eigen::Vector3d> maybe_normal )
{
    const Eigen::Vector3d normal = maybe_normal.value_or( triangleNormal( tri ) );
    const double ray_scaling = normal.dot( tri.v1 - ray.start_pos ) / normal.dot( ray.dir );
    LOG( LOG_TRACING ) << "| | Ray scaling: " << ray_scaling << std::endl;
    if( ray_scaling <= 0 ) return std::nullopt;
    const Eigen::Vector3d intersection_point =
        ray.start_pos + ray_scaling * ray.dir;

    const bool inside = normal.dot( ( tri.v2 - tri.v1 ).cross( intersection_point - tri.v1 ) ) >= 0 and
                        normal.dot( ( tri.v3 - tri.v2 ).cross( intersection_point - tri.v2 ) ) >= 0 and
                        normal.dot( ( tri.v1 - tri.v3 ).cross( intersection_point - tri.v3 ) ) >= 0;

    return inside ? std::optional<Eigen::Vector3d>( intersection_point ) : std::nullopt;
}

std::optional<TracePoint> traceRayOnTet( const cgogn::CMap3& map,
                                         const cgogn::CMap3::Volume& v,
                                         const Ray& ray,
                                         const std::vector<Normal>& normals )
{
    // The face on the input dart is the location that we start from.
    // Check all the other faces for intersection with the ray.
    LOG( LOG_TRACING ) << "Tracing on tet " << v << " from ray " << ray.start_pos.transpose() << " -> "
                       << ray.dir.transpose() << std::endl;
    bool found_it = false;
    std::pair<cgogn::CMap3::Face, Eigen::Vector3d> out;
    foreach_incident_edge( map, cgogn::CMap3::Face( v.dart_ ), [&]( cgogn::CMap3::Edge e ) {
        LOG( LOG_TRACING ) << "| Edge " << e << std::endl;
        const cgogn::CMap3::Face adj_face( phi2( map, e.dart_ ) );
        const Triangle tri = triangleOfFace( map, adj_face );
        LOG( LOG_TRACING ) << "| | Triangle: " << tri.v1.transpose() << " | " << tri.v2.transpose() << " | "
                           << tri.v3.transpose() << std::endl;

        const std::optional<Eigen::Vector3d> intersection =
            intersectionOf( ray, tri, normals.at( index_of( map, adj_face ) ).get( adj_face.dart_ ) );

        if( intersection.has_value() )
        {
            LOG( LOG_TRACING ) << "| | Found intersection: " << intersection.value().transpose() << std::endl;
            out.first = cgogn::CMap3::Face( phi3( map, adj_face.dart_ ) );
            out.second = intersection.value();
            found_it = true;
            return false;
        }
        return true;
    } );
    if( not found_it ) return std::nullopt;
    return out;
}

void tracingDebugOutput( const cgogn::CMap3& map,
                         const cgogn::CMap3::Volume& v,
                         const Ray& ray,
                         const SimplicialComplex& line,
                         SimplicialComplex& tets,
                         const size_t n )
{
    // Add the triangles of the current tet to tets
    foreach_incident_face( map, v, [&]( cgogn::CMap3::Face f ) {
        const Triangle from_face = triangleOfFace( map, f );
        const size_t next_vid = tets.points.size();
        tets.simplices.push_back( { next_vid + 0, next_vid + 1, next_vid + 2 } );
        tets.points.push_back( from_face.v1 );
        tets.points.push_back( from_face.v2 );
        tets.points.push_back( from_face.v3 );
        return true;
    } );

    // Create a simplicial complex with just one face: the triangle it's coming from
    const Triangle from_face = triangleOfFace( map, cgogn::CMap3::Face( v.dart_ ) );
    const SimplicialComplex from_face_complex( { { {0, 1, 2} }, { from_face.v1, from_face.v2, from_face.v3 } } );
    const io::VTKOutputObject from_face_output( from_face_complex );

    // Create a simplicial complex for the ray
    const SimplicialComplex ray_complex( { { {0} }, { ray.start_pos } } );
    io::VTKOutputObject ray_output( ray_complex );
    ray_output.addVertexField( "ray", ray.dir.transpose() );

    const io::VTKOutputObject line_output( line );
    const io::VTKOutputObject tets_output( tets );

    std::stringstream ss;
    ss << std::setw(3) << std::setfill('0') << n;
    std::string n_str(ss.str());

    io::outputSimplicialFieldToVTK( from_face_output, "from_face_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( ray_output, "ray_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( line_output, "line_" + n_str + ".vtu" );
    io::outputSimplicialFieldToVTK( tets_output, "tets_" + n_str + ".vtu" );
}

SimplicialComplex traceField( const cgogn::CMap3& map,
                              const cgogn::CMap3::Face& f,
                              const Eigen::Vector3d& start_point,
                              const Eigen::MatrixX3d& field,
                              const std::vector<Normal>& normals,
                              const bool debug_output )
{
    cgogn::CMap3::Face curr_face = f;
    Eigen::Vector3d curr_point = start_point;

    SimplicialComplex complex;
    SimplicialComplex debug_tets;
    size_t n = 0;

    complex.points.push_back( curr_point );
    while( not is_boundary( map, curr_face.dart_ ) )
    {
        const cgogn::CMap3::Volume curr_vol( curr_face.dart_ );
        const Ray search_ray( {curr_point, field.row( index_of( map, curr_vol ) )} );
        if( debug_output ) tracingDebugOutput( map, curr_vol, search_ray, complex, debug_tets, n++ );
        const std::optional<TracePoint> next_point = traceRayOnTet( map, curr_vol, search_ray, normals );
        if( not next_point.has_value() ) break;
        curr_face = next_point.value().first;
        curr_point = next_point.value().second;
        complex.points.push_back( curr_point );
        complex.simplices.push_back( { complex.points.size() - 2, complex.points.size() - 1 } );
    }

    return complex;
}