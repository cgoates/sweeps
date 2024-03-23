#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <QuadLayoutInput.hpp>
#include <Laplace.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <GlobalCellMarker.hpp>
#include <VTKOutput.hpp>
#include <Logging.hpp>
#include <Tracing.hpp>
#include <queue>

void foreachFaceWithVertsInSet( const topology::TetMeshCombinatorialMap& map,
                                const std::vector<bool>& set,
                                const std::function<bool( const topology::Face&, const size_t n )>& callback )
{
    const topology::CombinatorialMapBoundary bdry( map );
    const auto contains = [&]( const topology::Face& f ) {
        return set.at( map.vertexId( topology::Vertex( f.dart() ) ).id() ) and
               set.at( map.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) and
               set.at( map.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() );
    };
    const topology::CombinatorialMapRestriction base_surf( bdry, contains );
    size_t i = 0;
    iterateCellsWhile( base_surf, 2, [&]( const topology::Face& f ) { return callback( f, i++ ); } );
}

void foreachEdgeOnBoundaryOfSet( const topology::CombinatorialMapBoundary& map,
                                 const std::vector<bool>& set,
                                 const std::function<bool( const topology::Edge&, const size_t n )>& callback )
{
    const auto contains = [&]( const topology::Face& f ) {
        return set.at( map.vertexId( topology::Vertex( f.dart() ) ).id() ) and
               set.at( map.vertexId( topology::Vertex( phi( map, 1, f.dart() ).value() ) ).id() ) and
               set.at( map.vertexId( topology::Vertex( phi( map, -1, f.dart() ).value() ) ).id() );
    };
    const topology::CombinatorialMapRestriction base_surf( map, contains );
    const topology::CombinatorialMapBoundary edges( base_surf );
    size_t i = 0;
    iterateCellsWhile( edges, 1, [&]( const topology::Edge& e ) {
        const topology::Dart phi2 = phi( map, 2, e.dart() ).value();
        return callback( topology::Edge( phi2 ), i++ );
    } );
}

topology::Cell cellOfSimplex( const topology::TetMeshCombinatorialMap& map, const Simplex& s )
{
    const topology::Vertex v0 = map.vertexOfId( s.vertex( 0 ) );

    if( s.dim() == 0 ) return v0;

    SmallVector<VertexId, 4> vertices;
    for( size_t i = 0; i <= s.dim(); i++ ) vertices.push_back( s.vertex( i ) );

    const auto contains = [&]( const VertexId& v ) { return std::find( vertices.begin(), vertices.end(), v ) != vertices.end(); };

    std::optional<topology::Cell> out;
    iterateAdjacentCells( map, v0, s.dim(), [&]( const topology::Cell& c ) {
        const bool right_cell = iterateDartsOfCell( map, c, [&]( const topology::Dart& d ) {
            return contains( map.vertexId( topology::Vertex( d ) ) );
        } );
        if( right_cell ) out.emplace( c );
        return not right_cell;
    } );

    if( not out.has_value() ) throw std::runtime_error( "Bad simplex provided" );

    // FIXME: Somehow we need the coordinate permutation here as well
    if( s.dim() > 1 ) std::runtime_error( "Not implemented cellOfSimplex; needs coordinate permutation" );

    return out.value();
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

int main( int argc, char* argv[] )
{
    const std::vector<std::string> input_args(argv + 1, argv + argc);

    if( input_args.size() > 0 )
    {
        const SweepInput sweep_input =
            io::loadINPFile( SRC_HOME "/test/data/macaroni.inp", "Surface3", "Surface4" );

        if( input_args.front() == "transpose-barycoords" )
        {
            if( input_args.size() < 4 ) throw( "Bad input to transpose-barycoords" );

            io::rewriteBaryCoordFile( sweep_input.mesh, input_args.at( 1 ), input_args.at( 2 ), input_args.at( 3 ) );

            std::cout << "Successfully rewrote " << input_args.at( 2 ) << " to " << input_args.at( 3 ) << std::endl;
            return 0;
        }

        Timer t;

        std::cout << "output-laplace\n";
        t.start( 1 );
        const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
        std::cout << "tet mesh construct time: " << t.stop( 1 ) << std::endl;

        std::cout << "Built map\n";
        t.start( 3 );
        const std::vector<Normal> normals = faceNormals( map );
        std::cout << "Normals time: " << t.stop( 3 ) << std::endl;
        t.start( 4 );
        const Eigen::VectorXd ans =
            solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );
        std::cout << "Laplace time 2: " << t.stop( 4 ) << std::endl;

        const topology::CombinatorialMapBoundary bdry( map );

        const auto keep_face_sides = [&]( const topology::Face& f ) {
            return ( not sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) or
                        not sweep_input.zero_bcs.at(
                            bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) or
                        not sweep_input.zero_bcs.at(
                            bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() ) ) and
                    ( not sweep_input.one_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) or
                        not sweep_input.one_bcs.at(
                            bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) or
                        not sweep_input.one_bcs.at(
                            bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() ) );
        };

        const auto keep_face_base = [&]( const topology::Face& f ) {
            return sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( f.dart() ) ).id() ) and
                    sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ).id() ) and
                    sweep_input.zero_bcs.at( bdry.vertexId( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ).id() );
        };

        const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
        const topology::CombinatorialMapRestriction base( bdry, keep_face_base );

        const Eigen::Matrix3Xd grad = gradientsWithBoundaryCorrection( map, sides, ans, normals );

        const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ){
            return [&sweep_input, &map]( const topology::Vertex& v ) -> const Eigen::Vector3d& {
                return sweep_input.mesh.points.at( map.vertexId( v ).id() );
            };
        };

        if( std::find( input_args.begin(), input_args.end(), "output-laplace" ) != input_args.end() )
        {
            std::cout << "Outputting Laplace field...\n";
            io::VTKOutputObject output( sweep_input.mesh );
            output.addVertexField( "laplace", ans );
            output.addCellField( "gradients", grad.transpose() );
            io::outputSimplicialFieldToVTK( output, "test.vtu" );
        }

        if( input_args.front() == "inner-traces" )
        {
            t.start( 7 );
            std::cout << "Outputting centroid traces for all base faces...\n";
            SimplicialComplex all_lines2;
            foreachFaceWithVertsInSet(
                map, sweep_input.zero_bcs, [&]( const topology::Face& start_face, const size_t ) {
                    t.start( 9 );
                    const SimplicialComplex field_line =
                        traceField( map, start_face, centroid( map, start_face ), grad, normals );
                    append( all_lines2, field_line );
                    t.stop( 9 );
                    return true;
                } );

            io::VTKOutputObject output2( all_lines2 );
            io::outputSimplicialFieldToVTK( output2, "centroid_traces_01.vtu" );
            std::cout << "Time to trace each face: " << t.stop( 7 ) << " inner time: " << t.stop( 9 ) << std::endl;
        }
        if( input_args.front() == "outer-traces" )
        {
            std::cout << "Outputting all boundary traces...\n";

            SimplicialComplex boundary_lines;
            t.start( 7 );
            foreachEdgeOnBoundaryOfSet(
                bdry, sweep_input.zero_bcs, [&]( const topology::Edge& start_edge, const size_t ) {
                    t.start( 9 );
                    const SimplicialComplex field_line =
                        traceBoundaryField( sides, start_edge, 0.5, ans, vertex_positions( bdry ), false );
                    append( boundary_lines, field_line );
                    t.stop( 9 );
                    return true;
                } );
            std::cout << "Time to trace each edge: " << t.stop( 7 ) << " inner time: " << t.stop( 9 ) << std::endl;

            io::VTKOutputObject boundary_output( boundary_lines );
            io::outputSimplicialFieldToVTK( boundary_output, "boundary_lines_01.vtu" );
        }
        if( input_args.front() == "outer-surfaces" )
        {
            std::cout << "Generating outer surfaces\n";
            t.start( 6 );

            topology::GlobalCellMarker crossed_faces( sides, 2 );
            topology::GlobalCellMarker flooded_faces( sides, 2 );

            /*
                Mark every face that is crossed by a trace.
                Flood from every face that isn't marked, marking them and storing in a simplicial complex.
            */

            const auto is_marked = [&]( const topology::Face& f ) {
                return crossed_faces.isMarked( f ) or flooded_faces.isMarked( f );
            };
            const std::vector<std::vector<BarycentricPoint>> coords =
                io::loadBaryCoords( "/Users/caleb/Downloads/macaroni_layout_bary_coords_00" );

            SimplicialComplex boundary_lines;
            for( size_t i = 0; i < 40; i++ )
            {
                for( const auto& coord : { coords.at( i ).front(), coords.at( i ).back() } )
                {
                    const topology::Cell possible_trace_cell = cellOfSimplex( map, coord.simplex );
                    // FIXME: This shouldn't have to be an edge...
                    if( possible_trace_cell.dim() != 1 ) continue;
                    std::optional<topology::Cell> trace_cell;
                    iterateDartsOfCell( map, possible_trace_cell, [&]( const topology::Dart& d ) {
                        if( onBoundary( map, d ) and keep_face_sides( d ) )
                        {
                            trace_cell.emplace( topology::Cell( d, possible_trace_cell.dim() ) );
                            return false;
                        }
                        return true;
                    } );

                    if( not trace_cell.has_value() ) continue;

                    const bool edges_aligned = map.vertexId( possible_trace_cell.dart() ) == sides.vertexId( trace_cell.value().dart() );

                    const double b = edges_aligned ? coord.point( 1 ) : coord.point( 0 );

                    const SimplicialComplex field_line = traceBoundaryField(
                        sides, trace_cell.value(), b, ans, vertex_positions( bdry ), false, [&]( const topology::Face& f ) {
                            crossed_faces.mark( sides, f );
                        } );
                    //std::raise(SIGINT);
                    append( boundary_lines, field_line );
                }
            }

            io::VTKOutputObject boundary_output( boundary_lines );
            io::outputSimplicialFieldToVTK( boundary_output, "boundary_traces.vtu" );

            std::cout << "Done tracing\n";

            size_t k = 0;
            iterateCellsWhile( sides, 2, [&]( const topology::Face& f ) {
                if( not is_marked( f ) )
                {
                    // Flood from here
                    SimplicialComplex surface_component;
                    topology::flood2d(
                        sides,
                        f,
                        is_marked,
                        [&]( const topology::Face& f ) { flooded_faces.mark( sides, f ); },
                        [&]( const topology::Face& f ) {
                            // Convert face to triangle, add to simplicial complex, mark as visited
                            const auto tri = triangleOfFace( map, f );
                            addTriangleNoDuplicateChecking( surface_component, tri );
                        } );
                    io::VTKOutputObject surface_output( surface_component );
                    io::outputSimplicialFieldToVTK( surface_output, "outer_surface" + std::to_string( k++ ) + ".vtu" );
                }
                return true;
            } );
            std::cout << "Time for outside surfs: " << t.stop( 6 ) << std::endl;
        }
        else if( input_args.at( 0 ) == "inner-surfaces" )
        {
            const std::vector<std::vector<BarycentricPoint>> bary_curves =
                io::loadBaryCoords( "/Users/caleb/Downloads/macaroni_layout_bary_coords_00" );

            /*
                Mark every face that is crossed by a trace.
                Flood from every face that isn't marked, marking them and storing in a simplicial complex.
            */

            for( size_t i = 0; i < 40; i++ )
            {
                const std::vector<BarycentricPoint>& coords = bary_curves.at( i );
                std::cout << std::endl << i << std::endl;

                std::vector<std::vector<Eigen::Vector3d>> lines;
                SimplicialComplex sep_tris;

                size_t coord_ii = 0;
                for( const auto& coord : coords )
                {
                    bool boundary_of_base = false;
                    std::optional<topology::Cell> trace_cell;
                    iterateDartsOfCell( map, cellOfSimplex( map, coord.simplex ), [&]( const topology::Dart& d ) {
                        if( onBoundary( map, d ) )
                        {
                            const topology::Cell possible_cell( d, coord.simplex.dim() );
                            if( keep_face_sides( d ) )
                            {
                                boundary_of_base = true;
                                trace_cell.emplace( possible_cell );
                                return false;
                            }
                            else if( not ( coord.simplex.dim() < 2 and boundaryAdjacent( base, possible_cell ) ) )
                            {
                                trace_cell.emplace( possible_cell );
                                return false;
                            }
                        }
                        return true;
                    } );

                    if( not trace_cell.has_value() ) throw std::runtime_error( "No cell for barycentric coordinate!" );

                    const SimplicialComplex field_line = [&]() {
                        if( boundary_of_base )
                        {
                            const topology::Dart start_cell_dart =
                                trace_cell.value().dim() == 0 ? phi( bdry, 1, trace_cell.value().dart() ).value()
                                                              : trace_cell.value().dart();
                            const topology::Edge start_edge( start_cell_dart );

                            const double b = [&](){
                                if( trace_cell.value().dim() == 0 ) return 0.0;
                                else if( sides.vertexId( topology::Vertex( start_edge.dart() ) ) == coord.simplex.vertex( 0 ) )
                                    return coord.point( 1 );
                                else return coord.point( 0 );
                            }();

                            std::cout << "ON BOUNDARY " << i << "\n";

                            return traceBoundaryField( sides, start_edge, b, ans, vertex_positions( bdry ), false );
                        }
                        else
                        {
                            const Eigen::Vector3d pt =
                                expandBarycentric( map,
                                                   vertex_positions( map ),
                                                   trace_cell.value(),
                                                   coord );
                            return traceField( map, trace_cell.value(), pt, grad, normals, true );
                        }
                    }();

                    if( lines.size() > 0 )
                    {
                        if( field_line.points.size() > lines.back().size() )
                            accumulateTrianglesFromParallelLines( sep_tris, lines.back(), field_line.points );
                        else
                            accumulateTrianglesFromParallelLines( sep_tris, field_line.points, lines.back() );
                    }

                    if( boundary_of_base )
                    {
                        std::cout << "num_points: " << sep_tris.points.size() << " last point: " << field_line.points.back() << std::endl;
                        std::cout << "last tri: " << sep_tris.simplices.back() << std::endl;
                    }
                    lines.push_back( field_line.points );

                    coord_ii++;
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
