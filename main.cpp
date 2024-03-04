#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <QuadLayoutInput.hpp>
#include <Laplace.hpp>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/attributes.h>
#include <VTKOutput.hpp>
#include <Logging.hpp>
#include <Tracing.hpp>

void foreachFaceWithVertsInSet( const cgogn::CMap3& map,
                                const std::vector<bool>& set,
                                const std::function<bool( const cgogn::CMap3::Face&, const size_t n )>& callback )
{
    const auto contains = [&]( const VertexId& vid ) {
        return set.at( vid.id() );
    };
    size_t i = 0;
    foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
        const cgogn::Dart& d = f.dart_;
        const VertexId vid1 = index_of( map, cgogn::CMap3::Vertex( d ) );
        const VertexId vid2 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
        const VertexId vid3 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
        if( contains( vid1 ) and contains( vid2 ) and contains( vid3 ) ) return callback( f, i++ );
        return true;
    } );
}

std::optional<cgogn::CMap3::Edge> boundaryTracingEdge( const cgogn::CMap3& map,
                                                       const std::vector<bool>& set,
                                                       const cgogn::CMap3::Edge& start_e )
{
    const auto contains = [&]( const VertexId& vid ) { return set.at( vid.id() ); };

    std::optional<cgogn::CMap3::Edge> out = std::nullopt;
    foreach_dart_of_orbit( map, start_e, [&]( cgogn::Dart d ) {
        const cgogn::CMap3::Edge e( d );
        if( not is_boundary( map, d ) ) return true;
        const VertexId vid1 = index_of( map, cgogn::CMap3::Vertex( d ) );
        const VertexId vid2 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
        const VertexId vid3 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
        if( contains( vid1 ) and contains( vid2 ) and not contains( vid3 ) )
        {
            out.emplace( e );
            return false;
        }
        const VertexId vid4 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi<2, -1>( map, d ) ) );
        if( contains( vid1 ) and contains( vid2 ) and not contains( vid4 ) )
        {
            out.emplace( cgogn::CMap3::Edge( cgogn::phi2( map, d ) ) );
            return false;
        }
        return false;
    } );
    return out;
}

std::optional<cgogn::CMap3::Edge> boundaryTracingEdgeOnFace( const cgogn::CMap3& map,
                                                             const std::vector<bool>& set,
                                                             const cgogn::CMap3::Face& f )
{
    return boundaryTracingEdge( map, set, cgogn::CMap3::Edge( f.dart_ ) ).or_else( [&]() {
        return boundaryTracingEdge( map, set, cgogn::CMap3::Edge( phi1( map, f.dart_ ) ) ).or_else( [&]() {
            return boundaryTracingEdge( map, set, cgogn::CMap3::Edge( phi_1( map, f.dart_ ) ) );
        } );
    } );
}

void foreachEdgeOnBoundaryOfSet( const cgogn::CMap3& map,
                                 const std::vector<bool>& set,
                                 const std::function<bool( const cgogn::CMap3::Edge&, const size_t n )>& callback )
{
    size_t i = 0;
    foreach_cell( map, [&]( cgogn::CMap3::Edge e ) {
        const std::optional<cgogn::CMap3::Edge> boundary_e = boundaryTracingEdge( map, set, e );
        if( boundary_e.has_value() ) return callback( boundary_e.value(), i++ );
        return true;
    } );
}

void foreachBaryCoordOnSetBoundary( const cgogn::CMap3& map,
                                    const std::vector<bool>& set,
                                    std::vector<BarycentricPoint>& bary_coords,
                                    const std::function<bool( const cgogn::CMap3::Edge&, const double )>& callback )
{
    foreachEdgeOnBoundaryOfSet( map, set, [&]( const cgogn::CMap3::Edge& e, const size_t ) {
        const VertexId vid1 = index_of( map, cgogn::CMap3::Vertex( e.dart_ ) );
        const VertexId vid2 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, e.dart_ ) ) );
        const bool swap = vid2 < vid1;
        const Simplex s = swap ? Simplex( vid2, vid1 ) : Simplex( vid1, vid2 );
        std::cout << s << std::endl;
        const auto it = std::find_if(
            bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
        if( it != bary_coords.end() )
        {
            const auto pt = it->point;
            bary_coords.erase( it );
            return callback( e, swap ? pt( 0 ) : pt( 1 ) );
        }
        return true;
    } );
}

void foreachBaryCoordOnSet( const cgogn::CMap3& map,
                            const std::vector<bool>& set,
                            std::vector<BarycentricPoint>& bary_coords,
                            const std::function<bool( const cgogn::CMap3::Face&, const Eigen::Vector3d& )>& callback )
{
    foreachFaceWithVertsInSet( map, set, [&]( const cgogn::CMap3::Face& f, const size_t ) {
        /*
            1. Check if the face is in the barycoords. If so, call back on it.
            2. Check if each edge is in the barycoords.  If so, transfer to face and call back on it.
            3. Check if any vertex is in the barycoords.  If so, transfer to face and call back on it.
            4. Remove any of those barycoords from the list so that they aren't duplicated.
        */
        const cgogn::Dart& d = f.dart_;
        cgogn::Dart d1 = f.dart_;
        cgogn::Dart d2 = cgogn::phi1( map, d );
        cgogn::Dart d3 = cgogn::phi_1( map, d );
        VertexId vid1 = index_of( map, cgogn::CMap3::Vertex( d ) );
        VertexId vid2 = index_of( map, cgogn::CMap3::Vertex( d2 ) );
        VertexId vid3 = index_of( map, cgogn::CMap3::Vertex( d3 ) );

        const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
        Eigen::Vector3d pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
        Eigen::Vector3d pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d2 ) );
        Eigen::Vector3d pos3 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d3 ) );

        if( vid1 > vid2 )
        {
            std::swap( vid1, vid2 );
            std::swap( pos1, pos2 );
            std::swap( d1, d2 );
        }
        if( vid1 > vid3 )
        {
            std::swap( vid1, vid3 );
            std::swap( pos1, pos3 );
            std::swap( d1, d3 );
        }
        if( vid2 > vid3 )
        {
            std::swap( vid2, vid3 );
            std::swap( pos2, pos3 );
            std::swap( d2, d3 );
        }

        {
            const Simplex s( vid1, vid2, vid3 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                const Eigen::Vector3d pt = it->point( 0 ) * pos1 + it->point( 1 ) * pos2 + it->point( 2 ) * pos3;
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( f, pt ) ) return false;
            }
        }

        {
            const Simplex s( vid1 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( f, 0.99999 * pos1 + 5e-6 * pos2 + 5e-6 * pos3 ) ) return false;
            }
        }

        {
            const Simplex s( vid2 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( f, 5e-6 * pos1 + 0.99999 * pos2 + 5e-6 * pos3 ) ) return false;
            }
        }

        {
            const Simplex s( vid3 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( f, 5e-6 * pos1 + 5e-6 * pos2 + 0.99999 * pos3 ) ) return false;
            }
        }

        {
            const Simplex s( vid1, vid2 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::cout << "1,2|";
                const auto pt = it->point;
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( cgogn::CMap3::Face( d1 ), 0.999998 * pt(0) * pos1 + 0.999998 * pt(1) * pos2 + 2e-6 * pos3 ) ) return false;
            }
        }

        {
            const Simplex s( vid1, vid3 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::cout << "1,3|";
                const auto pt = it->point;
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( cgogn::CMap3::Face( d3 ), 0.999998 * pt(0) * pos1 + 0.999998 * pt(1) * pos3 + 2e-6 * pos2 ) ) return false;
            }
        }

        {
            const Simplex s( vid2, vid3 );
            const auto it = std::find_if(
                bary_coords.begin(), bary_coords.end(), [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
            if( it != bary_coords.end() )
            {
                std::cout << "2,3|";
                const auto pt = it->point;
                std::erase_if( bary_coords, [&]( const BarycentricPoint& pt ) { return pt.simplex == s; } );
                if( not callback( cgogn::CMap3::Face( d2 ), 0.999998 * pt(0) * pos2 + 0.999998 * pt(1) * pos3 + 2e-6 * pos1 ) ) return false;
            }
        }

        return true;
    } );
}

void append( SimplicialComplex& append_to, const SimplicialComplex& to_append )
{
    const size_t offset = append_to.points.size();
    append_to.points.insert( append_to.points.end(), to_append.points.begin(), to_append.points.end() );
    std::transform( to_append.simplices.cbegin(), to_append.simplices.cend(), std::back_inserter( append_to.simplices ), [&]( const Simplex& s ) {
        // FIXME: NOT GENERAL!!!
        return Simplex( s.vertex( 0 ).id() + offset, s.vertex( 1 ).id() + offset );
    } );
}

void accumulateTrianglesFromParallelLines( SimplicialComplex& append_to,
                                           const std::vector<Eigen::Vector3d>& short_line,
                                           const std::vector<Eigen::Vector3d>& long_line )
{
    const size_t short_offset = append_to.points.size();
    const size_t long_offset = short_offset + short_line.size();
    append_to.points.insert( append_to.points.end(), short_line.begin(), short_line.end() );
    append_to.points.insert( append_to.points.end(), long_line.begin(), long_line.end() );

    double size_ratio = (double)long_line.size() / (double)short_line.size();

    const auto short_point = [&]( const size_t i ) { return short_offset + i; };
    const auto long_point = [&]( const size_t j ) { return long_offset + j; };

    for( size_t j = 0, i = 1; i < short_line.size(); i++ )
    {
        append_to.simplices.push_back( Simplex( short_point( i - 1 ), short_point( i ), long_point( j++ ) ) );
        if( j < long_line.size() )
        {
            append_to.simplices.push_back( Simplex( short_point( i ), long_point( j - 1 ), long_point( j ) ) );
            if( j <= (size_t)std::ceil( i * size_ratio ) )
            {
                j++;
                if( j < long_line.size() ) append_to.simplices.push_back( Simplex( short_point( i ), long_point( j - 1 ), long_point( j ) ) );
            }
        }
    }
}

void flood2d( const cgogn::CMap3& map,
              const cgogn::CMap3::Face& f,
              const std::function<bool( const cgogn::CMap3::Face& )>& stop_condition,
              const std::function<void( const cgogn::CMap3::Face& )>& mark_callback,
              const std::function<void( const cgogn::CMap3::Face& )>& callback )
{
    std::queue<cgogn::CMap3::Face> to_process;
    to_process.push( f );

    for( ; not to_process.empty(); to_process.pop() )
    {
        const cgogn::CMap3::Face& curr_face = to_process.front();
        if( stop_condition( curr_face ) ) continue;
        callback( curr_face );
        mark_callback( curr_face );
        for( const auto& d :
             { phi2( map, curr_face.dart_ ), phi<1, 2>( map, curr_face.dart_ ), phi<-1, 2>( map, curr_face.dart_ ) } )
        {
            if( not stop_condition( cgogn::CMap3::Face( d ) ) )
            {
                to_process.push( cgogn::CMap3::Face( d ) );
            }
        }
    }
}

int main( int argc, char* argv[] )
{
    const std::vector<std::string> input_args(argv + 1, argv + argc);

    if( input_args.size() > 0 )
    {
        const SweepInput sweep_input =
            io::loadINPFile( "/Users/caleb/sweeps/attempt-sweep/test/data/macaroni.inp", "Surface3", "Surface4" );

        if( input_args.front() == "transpose-barycoords" )
        {
            if( input_args.size() < 4 ) throw( "Bad input to transpose-barycoords" );

            io::rewriteBaryCoordFile( sweep_input.mesh, input_args.at( 1 ), input_args.at( 2 ), input_args.at( 3 ) );

            std::cout << "Successfully rewrote " << input_args.at( 2 ) << " to " << input_args.at( 3 ) << std::endl;
            return 0;
        }

        cgogn::CMap3 map;
        mapFromInput( sweep_input.mesh, map );

        const std::vector<Normal> normals = faceNormals( map );
        const Eigen::VectorXd ans = solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

        const Eigen::MatrixX3d grad = gradients( map, ans, normals );

        if( std::find( input_args.begin(), input_args.end(), "output-laplace" ) != input_args.end() )
        {
            std::cout << "Outputting Laplace field...\n";
            io::VTKOutputObject output( sweep_input.mesh );
            output.addVertexField( "laplace", ans );
            output.addCellField( "gradients", grad );
            io::outputSimplicialFieldToVTK( output, "test.vtu" );
        }

        if( input_args.front() == "inner-traces" )
        {
            std::cout << "Outputting centroid traces for all base faces...\n";
            SimplicialComplex all_lines;
            foreachFaceWithVertsInSet( map, sweep_input.zero_bcs, [&]( const cgogn::CMap3::Face& start_face, const size_t ){
                const SimplicialComplex field_line = traceField( map, start_face, centroid( map, start_face ), grad, normals );
                append( all_lines, field_line );
                return true;
            } );

            io::VTKOutputObject output( all_lines );
            io::outputSimplicialFieldToVTK( output, "centroid_traces.vtu" );
        }
        if( input_args.front() == "outer-traces" )
        {
            std::cout << "Outputting all boundary traces...\n";
            SimplicialComplex boundary_lines;
            foreachEdgeOnBoundaryOfSet( map, sweep_input.zero_bcs, [&]( const cgogn::CMap3::Edge& start_edge, const size_t ) {
                const SimplicialComplex field_line = traceBoundaryField( map, start_edge, 0.5, ans, sweep_input.one_bcs, false );
                append( boundary_lines, field_line );
                return true;
            } );

            io::VTKOutputObject boundary_output( boundary_lines );
            io::outputSimplicialFieldToVTK( boundary_output, "boundary_lines.vtu" );
        }
        if( input_args.front() == "outer-surfaces" )
        {
            std::cout << "Generating outer surfaces\n";
            cgogn::CellMarker<cgogn::CMap3, cgogn::CMap3::Face> crossed_faces( map );
            cgogn::CellMarker<cgogn::CMap3, cgogn::CMap3::Face> source_or_target( map );
            cgogn::CellMarker<cgogn::CMap3, cgogn::CMap3::Face> flooded_faces( map );
            foreachFaceWithVertsInSet( map, sweep_input.zero_bcs, [&]( const cgogn::CMap3::Face& f, const auto ) {
                source_or_target.mark( f );
                return true;
            } );
            foreachFaceWithVertsInSet( map, sweep_input.one_bcs, [&]( const cgogn::CMap3::Face& f, const auto ) {
                source_or_target.mark( f );
                return true;
            } );

            /*
                Mark every face that is crossed by a trace.
                Flood from every face that isn't marked, marking them and storing in a simplicial complex.
            */

            const auto is_marked = [&]( cgogn::CMap3::Face f ) {
                return crossed_faces.is_marked( f ) or source_or_target.is_marked( f ) or flooded_faces.is_marked( f );
            };
            const std::vector<std::vector<BarycentricPoint>> coords =
                io::loadBaryCoords( "/Users/caleb/Downloads/macaroni_layout_bary_coords_00" );

            std::vector<BarycentricPoint> endpoints;
            for( size_t i = 0; i < 40; i++ )
            {
                endpoints.push_back( coords.at( i ).front() );
                endpoints.push_back( coords.at( i ).back() );
            }
            SimplicialComplex boundary_lines;
            foreachBaryCoordOnSetBoundary(
                map, sweep_input.zero_bcs, endpoints, [&]( const cgogn::CMap3::Edge& start_edge, const double b ) {
                    std::cout << "FOUND ONE EDGE\n";
                    const SimplicialComplex field_line = traceBoundaryField(
                        map, start_edge, b, ans, sweep_input.one_bcs, false, [&]( const cgogn::CMap3::Face& f ) {
                            crossed_faces.mark( f );
                        } );
                    append( boundary_lines, field_line );
                    return true;
                } );

            io::VTKOutputObject boundary_output( boundary_lines );
            io::outputSimplicialFieldToVTK( boundary_output, "boundary_traces.vtu" );

            size_t k = 0;
            foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
                f = cgogn::CMap3::Face( phi3( map, f.dart_ ) );
                if( is_boundary( map, f.dart_ ) and not is_marked( f ) )
                {
                    // Flood from here
                    SimplicialComplex surface_component;
                    flood2d(
                        map,
                        f,
                        is_marked,
                        [&]( const cgogn::CMap3::Face& f ) { flooded_faces.mark( f ); },
                        [&]( const cgogn::CMap3::Face& f ) {
                            // Convert face to triangle, add to simplicial complex, mark as visited
                            const auto tri = triangleOfFace( map, f );
                            addTriangleNoDuplicateChecking( surface_component, tri );
                        } );
                    io::VTKOutputObject surface_output( surface_component );
                    io::outputSimplicialFieldToVTK( surface_output, "outer_surface" + std::to_string( k++ ) + ".vtu" );
                }
                return true;
            } );
        }
        else if( input_args.at( 0 ) == "inner-surfaces" )
        {
            const std::vector<std::vector<BarycentricPoint>> bary_curves =
                io::loadBaryCoords( "/Users/caleb/Downloads/macaroni_layout_bary_coords_00" );
            for( size_t i = 0; i < 40; i++ )
            {
                const std::vector<BarycentricPoint>& coords = bary_curves.at( i );
                std::cout << std::endl << i << std::endl;

                std::vector<std::vector<Eigen::Vector3d>> lines;
                SimplicialComplex sep_tris;

                size_t coord_ii = 0;
                for( const auto& bc : coords )
                {
                    std::vector<BarycentricPoint> coord{ bc };
                    // FIXME: Iterate the points in the barycentric curve by looking for adjacencies
                    foreachBaryCoordOnSet( map, sweep_input.zero_bcs, coord, [&]( const cgogn::CMap3::Face& start_face, const Eigen::Vector3d& pt ) {

                        const double distance_from_major_axis = ( pt - Eigen::Vector3d( -15.5, 0, 0 ) ).norm();

                        const bool on_boundary = equals( distance_from_major_axis, 10.5, 1e-2 ) or
                            equals( distance_from_major_axis, 5, 1e-2 );


                        // const std::optional<cgogn::CMap3::Edge> boundary_e = boundaryTracingEdge( map, sweep_input.zero_bcs, cgogn::CMap3::Edge( start_face.dart_ ) );
                        const std::optional<cgogn::CMap3::Edge> boundary_e =
                            on_boundary
                            ? boundaryTracingEdgeOnFace( map, sweep_input.zero_bcs, start_face ) : std::nullopt;
                        const SimplicialComplex field_line = boundary_e.and_then( [&]( cgogn::CMap3::Edge start_edge ) -> std::optional<SimplicialComplex> {
                            const auto position =
                                cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
                            const Eigen::Vector3d pos1 =
                                cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( start_edge.dart_ ) );
                            const Eigen::Vector3d pos2 =
                                cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( phi1( map, start_edge.dart_ ) ) );

                            const double b = ( pt - pos1 ).norm() / ( pos2 - pos1 ).norm();

                            std::cout << "ON BOUNDARY " << i << "\n";

                            return traceBoundaryField( map, start_edge, b, ans, sweep_input.one_bcs, false );
                        } ).value_or( traceField( map, start_face, pt, grad, normals, true ) );

                        if( i != 1 or coord_ii != 2 )
                        {
                            if( lines.size() > 0 )
                            {
                                if( field_line.points.size() > lines.back().size() )
                                    accumulateTrianglesFromParallelLines( sep_tris, lines.back(), field_line.points );
                                else
                                    accumulateTrianglesFromParallelLines( sep_tris, field_line.points, lines.back() );
                            }

                            if( on_boundary )
                            {
                                std::cout << "num_points: " << sep_tris.points.size() << " last point: " << field_line.points.back() << std::endl;
                                std::cout << "last tri: " << sep_tris.simplices.back() << std::endl;
                            }
                            lines.push_back( field_line.points );
                        }

                        coord_ii++;
                        return true;
                    } );
                }

                io::VTKOutputObject boundary_output( sep_tris );
                io::outputSimplicialFieldToVTK( boundary_output, "hex_layout_skeleton" + std::to_string( i ) + ".vtu" );
            }
        }
    }
    else
    {
        std::cout << "No action specified" << std::endl;
    }
}
