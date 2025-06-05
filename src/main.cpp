#include <SimplexUtilities.hpp>
#include <MeshInput.hpp>
#include <QuadLayoutInput.hpp>
#include <Laplace.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <GlobalCellMarker.hpp>
#include <VTKOutput.hpp>
#include <OBJOutput.hpp>
#include <Logging.hpp>
#include <Tracing.hpp>
#include <queue>
#include <LevelSetCMap.hpp>
#include <DelaunayTriangulation.hpp>
#include <CriticalPoints.hpp>
#include "../test/SimplicialComplexTestCases.hpp"

using namespace reparam;

topology::Cell cellOfSimplex( const topology::TetMeshCombinatorialMap& map, const Simplex& s )
{
    const topology::Vertex v0 = map.vertexOfId( s.vertex( 0 ) );

    if( s.dim() == 0 ) return v0;

    SmallVector<VertexId, 4> vertices;
    for( size_t i = 0; i <= s.dim(); i++ ) vertices.push_back( s.vertex( i ) );

    const auto contains = [&]( const VertexId& v ) { return std::find( vertices.begin(), vertices.end(), v ) != vertices.end(); };

    const auto vertex_ids = indexingOrError( map, 0 );
    std::optional<topology::Cell> out;
    iterateAdjacentCells( map, v0, s.dim(), [&]( const topology::Cell& c ) {
        const bool right_cell = iterateDartsOfCell( map, c, [&]( const topology::Dart& d ) {
            return contains( VertexId( vertex_ids( topology::Vertex( d ) ) ) );
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
        switch( s.dim() )
        {
            case 0: return Simplex( s.vertex( 0 ).id() + offset );
            case 1: return Simplex( s.vertex( 0 ).id() + offset, s.vertex( 1 ).id() + offset );
            case 2:
                return Simplex( s.vertex( 0 ).id() + offset, s.vertex( 1 ).id() + offset, s.vertex( 2 ).id() + offset );
            case 3:
                return Simplex( s.vertex( 0 ).id() + offset,
                                s.vertex( 1 ).id() + offset,
                                s.vertex( 2 ).id() + offset,
                                s.vertex( 3 ).id() + offset );
            default: throw std::runtime_error( "Unsupported dimension" );
        }
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

// Check Delaunay criterion for a single tetrahedron
bool isDelaunayTetrahedron( const Eigen::MatrixXd& vertices, const Eigen::Vector4i& tet )
{
    // Extract tetrahedron vertices
    const auto tet_vertices = vertices( tet, Eigen::all );

    const Tetrahedron tetra{ tet_vertices.row( 0 ), tet_vertices.row( 1 ), tet_vertices.row( 2 ), tet_vertices.row( 3 ) };

    // Compute circumcenter
    const Eigen::Vector3d circumc = circumcenter( tetra );
    const double circumradius = ( circumc - tet_vertices.row( 0 ).transpose() ).norm();

    // Check if any other point is inside circumsphere
    for( int i = 0; i < vertices.rows(); ++i )
    {
        // Skip vertices of the current tetrahedron
        if( std::find( tet.data(), tet.data() + 4, i ) != tet.data() + 4 )
        {
            continue;
        }

        // Check if point is inside circumsphere
        if( ( vertices.row( i ).transpose() - circumc ).norm() < circumradius )
        {
            return false;
        }
    }

    return true;
}

// Verify Delaunay criterion for all tetrahedra
std::vector<int> verifyDelaunayTetrahedralization( const SimplicialComplex& mesh )
{
    // Convert output to Eigen matrices for easier manipulation
    Eigen::MatrixXd vertices( mesh.points.size(), 3 );
    for( int i = 0; i < vertices.rows(); ++i )
    {
        vertices.row( i ) = mesh.points.at( i ).transpose();
    }

    std::vector<int> failed_tetrahedra;

    // Check each tetrahedron
    for( size_t i = 0; i < mesh.simplices.size(); ++i )
    {
        Eigen::Vector4i tet( mesh.simplices.at( i ).vertex( 0 ).id(),
                             mesh.simplices.at( i ).vertex( 1 ).id(),
                             mesh.simplices.at( i ).vertex( 2 ).id(),
                             mesh.simplices.at( i ).vertex( 3 ).id() );

        if( !isDelaunayTetrahedron( vertices, tet ) )
        {
            failed_tetrahedra.push_back( i );
        }
    }

    return failed_tetrahedra;
}

int main( int argc, char* argv[] )
{
    const std::vector<std::string> input_args(argv + 1, argv + argc);

    if( input_args.size() > 0 )
    {
        const SweepInput sweep_input = [&]() {
            if( std::find( input_args.begin(), input_args.end(), "femur" ) != input_args.end() )
                return SweepInputTestCases::femur();
            else if( std::find( input_args.begin(), input_args.end(), "hook" ) != input_args.end() )
                return SweepInputTestCases::hook();
            else if( std::find( input_args.begin(), input_args.end(), "ventricle" ) != input_args.end() )
                return SweepInputTestCases::ventricle();
            else if( std::find( input_args.begin(), input_args.end(), "flange" ) != input_args.end() )
                return SweepInputTestCases::flange();
            else if( std::find( input_args.begin(), input_args.end(), "cube" ) != input_args.end() )
                return SweepInputTestCases::twelveTetCube();
            else if( std::find( input_args.begin(), input_args.end(), "bunny" ) != input_args.end() )
                return SweepInputTestCases::bunny();
            else if( std::find( input_args.begin(), input_args.end(), "spring" ) != input_args.end() )
                return SweepInputTestCases::spring();
            else if( std::find( input_args.begin(), input_args.end(), "half_bullet" ) != input_args.end() )
                return SweepInputTestCases::bullet();
            else if( std::find( input_args.begin(), input_args.end(), "bullet" ) != input_args.end() )
                return SweepInputTestCases::bullet_full();
            else if( std::find( input_args.begin(), input_args.end(), "part_torus_in_sphere" ) != input_args.end() )
                return SweepInputTestCases::part_torus_in_sphere();
            else if( std::find( input_args.begin(), input_args.end(), "macaroni_coarse" ) != input_args.end() )
                return io::loadINPFile( SRC_HOME "/test/data/macaroni_coarse.inp", "Surface3", "Surface4" );
            else if( std::find( input_args.begin(), input_args.end(), "pot_counter" ) != input_args.end() )
                return SweepInputTestCases::pot_counter();
            else if( std::find( input_args.begin(), input_args.end(), "pot_counter_remeshed" ) != input_args.end() )
                return SweepInputTestCases::pot_counter( true );
            else if( std::find( input_args.begin(), input_args.end(), "sphere" ) != input_args.end() )
                return SweepInputTestCases::neighborhood();
            else if( std::find( input_args.begin(), input_args.end(), "counter2" ) != input_args.end() )
                return io::loadINPFile( SRC_HOME "/test/data/counter2.inp", "Surface1", "Surface6" );
            else
                return SweepInputTestCases::macaroni();
        }();

        if( input_args.front() == "transpose-barycoords" )
        {
            if( input_args.size() < 4 ) throw std::runtime_error( "Bad input to transpose-barycoords" );

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

        const topology::CombinatorialMapBoundary bdry( map );

        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
        const auto vertex_ids = indexingOrError( map, 0 );

        const auto keep_face_sides = [&]( const topology::Face& f ) {
            return ( not sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) or
                        not sweep_input.zero_bcs.at(
                            bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) or
                        not sweep_input.zero_bcs.at(
                            bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) ) ) and
                    ( not sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) or
                        not sweep_input.one_bcs.at(
                            bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) or
                        not sweep_input.one_bcs.at(
                            bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) ) );
        };

        const auto keep_face_base = [&]( const topology::Face& f ) {
            return sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
                    sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
                    sweep_input.zero_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
        };

        const auto keep_face_target = [&]( const topology::Face& f ) {
            return sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( f.dart() ) ) ) and
                   sweep_input.one_bcs.at(
                       bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
                   sweep_input.one_bcs.at(
                       bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
        };

        const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
        const topology::CombinatorialMapRestriction base( bdry, keep_face_base );
        const topology::CombinatorialMapRestriction target( bdry, keep_face_target );

        const auto sides_vertex_ids = indexingOrError( sides, 0 );
        const auto base_vertex_ids = indexingOrError( base, 0 );

        const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ){
            const auto vertex_ids = indexingOrError( map, 0 );
            return [&sweep_input, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
                return sweep_input.mesh.points.at( vertex_ids( v ) );
            };
        };

        if( input_args.at( 0 ) == "output-setup" )
        {
            {
                SimplicialComplex one_level_set;
                addAllTriangles( one_level_set, base, vertex_positions( base ) );
                io::VTKOutputObject output( one_level_set );
                io::outputSimplicialFieldToVTK( output, "base.vtu" );
            }
            {
                SimplicialComplex one_level_set;
                addAllTriangles( one_level_set, target, vertex_positions( target ) );
                io::VTKOutputObject output( one_level_set );
                io::outputSimplicialFieldToVTK( output, "target.vtu" );
            }
            io::VTKOutputObject output( sweep_input.mesh );
            output.addVertexField( "vertex_ids", Eigen::VectorXd::LinSpaced( sweep_input.mesh.points.size(), 0, sweep_input.mesh.points.size() - 1 ) );
            io::outputSimplicialFieldToVTK( output, "mesh.vtu" );

            if( std::find( input_args.begin(), input_args.end(), "edge-weights" ) != input_args.end() )
            {
                size_t n_interior = 0;
                size_t n_boundary = 0;
                const std::vector<double> edge_weights = edgeWeightsLaplace3d(
                    map, vertex_positions( map ), faceNormals( map ), Laplace3dEdgeWeights::BarycentricDual );
                const auto edge_ids = indexingOrError( map, 1 );
                std::vector<topology::Edge> negative_weight_edges;
                iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
                    const size_t eid = edge_ids( e );
                    if( edge_weights.at( eid ) < 0 )
                    {
                        negative_weight_edges.push_back( e );
                        std::cout << ( boundaryAdjacent( map, e ) ? "Boundary edge with negative weight "
                                                                  : "Interior edge with negative weight " )
                                  << e << std::endl;
                        if( boundaryAdjacent( map, e ) )
                            n_boundary++;
                        else
                            n_interior++;

                        io::outputDualFace(
                            map, vertex_positions( map ), e, std::format( "{:04d}", negative_weight_edges.size() ) );
                    }
                    return true;
                } );
                std::cout << n_boundary << " boundary edges with negative weight\n";
                std::cout << n_interior << " interior edges with negative weight\n";
                io::outputEdgeChain( map, vertex_positions( map ), negative_weight_edges, "negative_weight_edges.vtu" );
            }

            if( std::find( input_args.begin(), input_args.end(), "delaunay" ) != input_args.end() )
            {
                std::vector<int> failed_tetrahedra = verifyDelaunayTetrahedralization( sweep_input.mesh );
                if( failed_tetrahedra.empty() )
                {
                    std::cout << "Delaunay criterion satisfied for all tetrahedra\n";
                }
                else
                {
                    std::cout << "Delaunay criterion violated for " << failed_tetrahedra.size() << " tetrahedra\n";
                }
            }

            double edge_length = 0;
            double min_edge_length = std::numeric_limits<double>::infinity();
            double max_edge_length = -std::numeric_limits<double>::infinity();
            iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
                const double length = edgeLength( map, vertex_positions( map ), e );
                edge_length += length;
                min_edge_length = std::min( min_edge_length, length );
                max_edge_length = std::max( max_edge_length, length );
                return true;
            } );
            edge_length /= cellCount( map, 1 );

            std::cout << "Average edge length: " << edge_length << std::endl;
            std::cout << "Min edge length: " << min_edge_length << std::endl;
            std::cout << "Max edge length: " << max_edge_length << std::endl;

            return 0;
        }

        t.start( 3 );
        const std::vector<Normal> normals = faceNormals( map );
        std::cout << "Normals time: " << t.stop( 3 ) << std::endl;
        t.start( 4 );
        const Eigen::VectorXd ans =
            sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );
        std::cout << "Laplace time 2: " << t.stop( 4 ) << std::endl;

        const Eigen::Matrix3Xd grad = gradientsWithBoundaryCorrection( map, sides, ans, normals );

        if( std::find( input_args.begin(), input_args.end(), "output-laplace" ) != input_args.end() )
        {
            std::cout << "Outputting Laplace field...\n";
            io::VTKOutputObject output( sweep_input.mesh );
            output.addVertexField( "laplace", ans );
            output.addCellField( "gradients", grad.transpose() );
            io::outputSimplicialFieldToVTK( output, "test.vtu" );
        }

        if( std::find( input_args.begin(), input_args.end(), "output-critical-points" ) != input_args.end() )
        {
            std::cout << "Outputting critical points...\n";
            SimplicialComplex critical_points;
            const auto func = [&]( const topology::Vertex& v ) {
                return ans( vertex_ids( v ) );
            };
            std::vector<int> euler_characteristics;
            std::vector<double> laplace_values;
            iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
                if( sweep_input.zero_bcs.at( vertex_ids( v ) ) ) return true;
                const int euler = reparam::lowerLinkEulerCharacteristic( map, v, func );
                if( euler != 1 )
                {
                    critical_points.points.push_back( vertex_positions( map )( v ) );
                    critical_points.simplices.push_back( Simplex( critical_points.points.size() - 1 ) );
                    euler_characteristics.push_back( euler );
                    laplace_values.push_back( ans( vertex_ids( v ) ) );
                }
                return true;
            } );
            io::VTKOutputObject output( critical_points );
            const auto euler_field = Eigen::Map<Eigen::VectorXi>( euler_characteristics.data(), euler_characteristics.size(), 1 );
            output.addVertexField( "euler", euler_field.cast<double>() );
            output.addVertexField( "laplace", Eigen::Map<Eigen::VectorXd>( laplace_values.data(), laplace_values.size(), 1 ) );
            io::outputSimplicialFieldToVTK( output, "critical_points.vtu" );
        }

        if( input_args.front() == "inner-traces" )
        {
            t.start( 7 );
            std::cout << "Outputting centroid traces for all base faces...\n";
            SimplicialComplex all_lines2;
            iterateCellsWhile( base, 2, [&]( const topology::Face& f ) {
                const topology::Face start_face = bdry.toUnderlyingCell( f );
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
            iterateCellsWhile( base, 1, [&]( const topology::Edge& e ) {
                if( not boundaryAdjacent( base, e ) ) return true;
                const topology::Edge start_edge( phi( bdry, 2, e.dart() ).value() );
                t.start( 9 );
                const reparam::Trace field_line =
                    traceBoundaryField( sides, start_edge, 0.5, ans, vertex_positions( bdry ), false );
                append( boundary_lines, field_line.mComplex );
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
                io::loadBaryCoords( SRC_HOME "/test/data/macaroni_layout_bary_coords_00" );

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

                    const bool edges_aligned = vertex_ids( topology::Vertex( possible_trace_cell.dart() ) ) ==
                                               sides_vertex_ids( topology::Vertex( trace_cell.value().dart() ) );

                    const double b = edges_aligned ? coord.point( 1 ) : coord.point( 0 );

                    const reparam::Trace field_line =
                        traceBoundaryField( sides, trace_cell.value(), b, ans, vertex_positions( bdry ), false );
                    for( const topology::Cell& f : field_line.mBaseCells )
                    {
                        crossed_faces.mark( sides, f );
                    }
                    append( boundary_lines, field_line.mComplex );
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
                io::loadBaryCoords( SRC_HOME "/test/data/macaroni_layout_bary_coords_00" );

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
                                else if( sides_vertex_ids( topology::Vertex( start_edge.dart() ) ) == coord.simplex.vertex( 0 ).id() )
                                    return coord.point( 1 );
                                else return coord.point( 0 );
                            }();

                            std::cout << "ON BOUNDARY " << i << "\n";

                            return traceBoundaryField( sides, start_edge, b, ans, vertex_positions( bdry ), false )
                                .mComplex;
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
                }

                io::VTKOutputObject boundary_output( sep_tris );
                io::outputSimplicialFieldToVTK( boundary_output, "hex_layout_skeleton" + std::to_string( i ) + ".vtu" );
            }
        }
        else if( input_args.at( 0 ) == "levelsets" )
        {
            SimplicialComplex all_level_sets;

            const size_t n_level_sets = 30;
            Eigen::VectorXd values = Eigen::VectorXd::LinSpaced( n_level_sets, 0.0, 1.0 );

            const auto func = [&]( const topology::Vertex& v ) {
                return ans( vertex_ids( v ) );
            };

            const bool output_individual_level_sets = std::find( input_args.begin(), input_args.end(), "individual" ) != input_args.end();

            addAllTriangles( all_level_sets, base, vertex_positions( base ) );
            if( output_individual_level_sets )
            {
                SimplicialComplex one_level_set;
                addAllTriangles( one_level_set, base, vertex_positions( base ) );
                io::VTKOutputObject output( one_level_set );
                io::outputSimplicialFieldToVTK( output, "level_set_0.vtu" );
            }
            for( Eigen::Index i = 1; i < values.size() - 1; i++ )
            {
                const double val = values( i );
                const auto level = std::make_shared<const topology::LevelSetCMap>( map, func, val );
                const auto level_positions = levelSetVertexPositions( *level, vertex_positions( map ) );
                const topology::DelaunayTriangulation level_tri( level, level_positions );
                const auto level_tri_positions = delaunayTriangulationVertexPositions( level_tri, level_positions );

                addAllTriangles( all_level_sets, level_tri, level_tri_positions );

                if( output_individual_level_sets )
                {
                    SimplicialComplex one_level_set;
                    addAllTriangles( one_level_set, level_tri, level_tri_positions );
                    io::VTKOutputObject output( one_level_set );
                    io::outputSimplicialFieldToVTK( output, "level_set_" + std::to_string( i ) + ".vtu" );
                }
            }
            addAllTriangles( all_level_sets, target, vertex_positions( target ) );
            if( output_individual_level_sets )
            {
                SimplicialComplex one_level_set;
                addAllTriangles( one_level_set, target, vertex_positions( target ) );
                io::VTKOutputObject output( one_level_set );
                io::outputSimplicialFieldToVTK( output, "level_set_" + std::to_string( values.size() - 1 ) + ".vtu" );
            }

            io::VTKOutputObject output( all_level_sets );
            io::outputSimplicialFieldToVTK( output, "level_sets.vtu" );
        }
        else if( input_args.at( 0 ) == "parameterize" )
        {
            /// LOAD 2D PARAMETERIZATION

            const auto bdry_positions = vertex_positions( bdry );
            size_t n_errors = 0;
            SimplicialComplex error_verts;
            // const auto negative_field = -1.0 * ans;
            topology::GlobalCellMarker traced_vertices( map, 0 );
            Eigen::Matrix3Xd param = Eigen::MatrixXd::Zero( 3, cellCount( map, 0 ) );

            std::cout << "Parameterize base\n";
            const auto base_pos = vertex_positions( base );
            iterateCellsWhile( base, 0, [&]( const topology::Vertex& v ) {
                traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
                param.col( base_vertex_ids( v ) ) = base_pos( v );
                return true;
            } );

            const auto reverse_ans = -1 * ans;
            const auto reverse_grad = -1 * grad;

            std::cout << "Parameterize sides\n";
            iterateCellsWhile( sides, 0, [&]( const topology::Vertex& v ) {
                if( traced_vertices.isMarked( bdry.toUnderlyingCell( v ) ) ) return true;
                traced_vertices.mark( map, bdry.toUnderlyingCell( v ) );
                try{
                    const SimplicialComplex line =
                        traceBoundaryField( sides, v, 1.0, reverse_ans, bdry_positions, false ).mComplex;
                    param.col( sides_vertex_ids( v ) ).head( 2 ) = line.points.back().head( 2 );
                    param( 2, sides_vertex_ids( v ) ) = ans( sides_vertex_ids( v ) );
                }
                catch( const std::runtime_error& e )
                {
                    std::cerr << "Skipping " << vertex_ids( v ) << " " << v << " with exception " << e.what();
                    n_errors++;
                    error_verts.points.push_back( bdry_positions( v ) );
                    error_verts.simplices.push_back( Simplex( error_verts.points.size() - 1 ) );
                    std::cout << " " << bdry_positions( v ).transpose() << std::endl;
                    throw e;
                }
                catch( const std::out_of_range& e )
                {
                    std::cerr << "Skipping " << vertex_ids( v ) << " " << v << " with exception " << e.what();
                    n_errors++;
                    error_verts.points.push_back( bdry_positions( v ) );
                    error_verts.simplices.push_back( Simplex( error_verts.points.size() - 1 ) );
                    std::cout << " " << bdry_positions( v ).transpose() << std::endl;
                }
                return true;
            } );

            std::cout << "Parameterize interior\n";
            const auto pos = vertex_positions( map );
            iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
                if( traced_vertices.isMarked( v ) ) return true;
                try
                {
                    const SimplicialComplex line = traceField( map, v, pos( v ), reverse_grad, normals, false );
                    param.col( vertex_ids( v ) ).head( 2 ) = line.points.back().head( 2 );
                    param( 2, vertex_ids( v ) ) = ans( vertex_ids( v ) );
                }
                catch( const std::runtime_error& e )
                {
                    std::cerr << "Skipping " << vertex_ids( v ) << " " << v << " with exception " << e.what();
                    n_errors++;
                    error_verts.points.push_back( pos( v ) );
                    error_verts.simplices.push_back( Simplex( error_verts.points.size() - 1 ) );
                    std::cout << " " << pos( v ).transpose() << std::endl;
                    throw e;
                }
                catch( const std::out_of_range& e )
                {
                    std::cerr << "Skipping " << vertex_ids( v ) << " " << v << " with exception " << e.what();
                    n_errors++;
                    error_verts.points.push_back( pos( v ) );
                    error_verts.simplices.push_back( Simplex( error_verts.points.size() - 1 ) );
                    std::cout << " " << pos( v ).transpose() << std::endl;
                }
                return true;
            } );

            io::VTKOutputObject output2( error_verts );
            io::outputSimplicialFieldToVTK( output2, "error_verts.vtu" );
            const auto param_positions = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
                return param.col( vertex_ids( v ) );
            };
            SimplicialComplex inverted_tets;
            SimplicialComplex interted_traces;
            iterateCellsWhile( map, 3, [&]( const topology::Volume& vol ) {
                if( isInverted( map, vol, param_positions ) )
                {
                    std::cerr << "Inverted volume " << vol << std::endl;
                    addTetNoDuplicateChecking( inverted_tets, map, pos, vol );
                    iterateAdjacentCells( map, vol, 0, [&]( const topology::Vertex& v ) {
                        const bool on_bdry = boundaryAdjacent( map, v );
                        std::cerr << "On boundary? " << on_bdry << std::endl;

                        std::cout << "Tracing from vertex of inverted tet...\n";
                        if( on_bdry )
                        {
                            // get the equvalent boundary vertex
                            iterateDartsOfCell( map, v, [&]( const topology::Dart& d ) {
                                if( onBoundary( map, d ) )
                                {
                                    const topology::Vertex bdry_v( phi( sides, 1, d ).value() );
                                    const reparam::Trace line =
                                        traceBoundaryField( sides, bdry_v, 1.0, reverse_ans, bdry_positions, false );
                                    append( interted_traces, line.mComplex );
                                    return false;
                                }
                                return true;
                            } );
                        }
                        else
                        {
                            const SimplicialComplex line = traceField( map, v, pos( v ), reverse_grad, normals, false );
                            append( interted_traces, line );
                        }

                        return true;
                    } );
                }
                return true;
            } );

            io::VTKOutputObject output_inverted_tets( inverted_tets );
            io::outputSimplicialFieldToVTK( output_inverted_tets, "inverted_tets.vtu" );

            io::VTKOutputObject output_inverted_traces( interted_traces );
            io::outputSimplicialFieldToVTK( output_inverted_traces, "inverted_traces.vtu" );

            std::cout << "Outputting Parameterization...\n";
            std::cout << n_errors << " errors in tracing\n";
            io::VTKOutputObject output( sweep_input.mesh );
            output.addVertexField( "param", param.transpose() );
            output.addVertexField( "laplace", ans );
            io::outputSimplicialFieldToVTK( output, "param.vtu" );
        }
        else if( input_args.at( 0 ) == "output-bdry-obj" )
        {
            io::outputTetMeshBoundaryToOBJ( map, "bdry.obj" );
        }
    }
    else
    {
        std::cout << "No action specified" << std::endl;
    }
}
