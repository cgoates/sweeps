#include <Foliation.hpp>
#include <Tracing.hpp>
#include <Cell.cpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <LevelSetCMap.hpp>
#include <Logging.hpp>
#include <optional>

#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <DelaunayTriangulation.hpp>
#include <ReversedCombinatorialMap.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshCircleMapping.hpp>
#include <TriangleMeshMapping.hpp>
#include <Laplace.hpp>
#include <CommonUtils.hpp>
#include <SplineSpace.hpp>
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>

namespace reparam
{
    std::vector<TraceLevelSetIntersection> levelSetIntersections( const Trace& trace,
                                                                  const topology::CombinatorialMap& traced_map,
                                                                  const std::vector<double>& level_set_values )
    {
        std::vector<TraceLevelSetIntersection> out;
        out.reserve( level_set_values.size() );

        const auto& vert_values = trace.mHarmonicFuncValues;

        size_t trace_vert_ii = 1;
        for( const double value : level_set_values )
        {
            while( vert_values.at( trace_vert_ii ) < value ) trace_vert_ii++;

            const auto [s, t] =
                inverseLinear( vert_values.at( trace_vert_ii - 1 ), vert_values.at( trace_vert_ii ), value );

            const Eigen::Vector3d pt =
                s * trace.mComplex.points.at( trace_vert_ii - 1 ) + t * trace.mComplex.points.at( trace_vert_ii );

            const topology::Cell& base_cell = trace.mBaseCells.at( trace_vert_ii - 1 );

            out.push_back( { pt, indexingOrError( traced_map, base_cell.dim() )( base_cell ) } );
        }

        return out;
    }

    std::map<topology::Vertex, double>
        thetaValues( const topology::CombinatorialMap& level_set,
                     const VertexPositionsFunc& level_set_positions,
                     const std::function<size_t( const topology::Edge& )>& underlying_face_id_of_edge,
                     const TraceLevelSetIntersection& intersection )
    {
        const topology::CombinatorialMapBoundary bdry( level_set );
        const auto bdry_positions = boundaryVertexPositions( bdry, level_set_positions );

        const topology::Dart start_d = [&bdry]() {
            topology::Dart d;
            iterateDartsWhile( bdry, [&]( const topology::Dart& a ) {
                d = a;
                return false;
            } );
            return d;
        }();

        const auto is_intersection_edge = [&]( const topology::Edge& bdry_e ) {
            const topology::Edge level_e = bdry.toUnderlyingCell( bdry_e );
            return ( underlying_face_id_of_edge( level_e ) == intersection.second );
        };

        std::map<topology::Vertex, double> out;

        topology::Dart d = start_d;
        double cumulative_length = 0.0;
        std::optional<double> intersection_point = std::nullopt;
        do
        {
            out.insert( { bdry.toUnderlyingCell( topology::Vertex( d ) ), cumulative_length } );

            if( not intersection_point and is_intersection_edge( d ) )
            {
                intersection_point.emplace( cumulative_length + ( intersection.first - bdry_positions( d ) ).norm() );
            }

            cumulative_length += edgeLength( bdry, bdry_positions, d );
            d = phi( bdry, 1, d ).value();
        } while( d != start_d );

        if( not intersection_point ) throw std::runtime_error( "Trace doesn't interact with the level set!" );

        const double factor = 2 * std::numbers::pi / cumulative_length;
        for( auto& pr : out )
        {
            pr.second = 2 * std::numbers::pi - ( pr.second - intersection_point.value() ) * factor;
        }

        return out;
    }


    std::map<topology::Vertex, Eigen::Vector2d>
        boundaryConstraints( const topology::CombinatorialMap& cut_cmap,
                             const VertexPositionsFunc& cut_cmap_positions,
                             const size_t n_cuts,
                             const std::function<bool( const topology::Vertex& )>& is_cut_extremity,
                             const std::function<bool( const topology::Vertex& )>& is_start_v )
    {
        if( cut_cmap.dim() != 2 ) throw std::invalid_argument( "boundaryConstraints only supports 2d cmaps" );
        const topology::CombinatorialMapBoundary bdry( cut_cmap );
        const auto bdry_positions = boundaryVertexPositions( bdry, cut_cmap_positions );
        const auto bdry_vert_ids = indexingOrError( bdry, 0 );

        const topology::Dart start_d = [&]() {
            std::optional<topology::Dart> d;
            iterateDartsWhile( bdry, [&]( const topology::Dart& a ) {
                if( is_start_v( bdry.toUnderlyingCell( topology::Vertex( a ) ) ) )
                {
                    d.emplace( a );
                    return false;
                }
                return true;
            } );
            if( not d.has_value() ) throw std::invalid_argument( "No vertex found for which is_start_v is true" );
            return d.value();
        }();

        std::map<topology::Vertex, Eigen::Vector2d> out;

        topology::Dart d = start_d;

        const size_t n_sides = n_cuts * 4;

        const std::vector<Eigen::Vector2d> ngon_verts = util::regularNGonVertices( n_sides );

        for( size_t side_ii = 0; side_ii < n_sides; side_ii++ )
        {
            std::map<topology::Vertex, double> side_positions;
            double cumulative_length = 0.0;
            do
            {
                side_positions.insert( { bdry.toUnderlyingCell( topology::Vertex( d ) ), cumulative_length } );
                cumulative_length += edgeLength( bdry, bdry_positions, d );
                d = phi( bdry, 1, d ).value();
            } while( not is_cut_extremity( bdry.toUnderlyingCell( topology::Vertex( d ) ) ) );

            const double factor = 1 / cumulative_length;
            for( auto& pr : side_positions )
            {
                const double s = pr.second * factor;
                out.emplace( pr.first, ( 1.0 - s ) * ngon_verts.at( side_ii ) + s * ngon_verts.at( side_ii + 1 ) );
            }
        }

        return out;
    }

    constexpr bool log_level_set_based_tracing = false;
    void levelSetBasedTracing( const SweepInput& sweep_input,
                               const std::vector<double> level_set_values,
                               const std::function<void( const std::vector<FoliationLeaf>& )>& callback )
    {
        const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
        const std::vector<Normal> normals = faceNormals( map );
        const Eigen::VectorXd sol = reparam::sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

        for( size_t ii = 1; ii < level_set_values.size() - 1; ii++ )
        {
            const double val = level_set_values.at( ii );
            if( ( sol.array() == val ).any() )
                throw std::runtime_error( "Cannot have a level set directly on a vertex" );
        }

        LOG( log_level_set_based_tracing ) << "FINISHED LAPLACE\n\n";

        const topology::CombinatorialMapBoundary bdry( map );

        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );
        const auto map_face_ids = indexingOrError( map, 2 );

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
                sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, 1, f.dart() ).value() ) ) ) and
                sweep_input.one_bcs.at( bdry_vertex_ids( topology::Vertex( phi( bdry, -1, f.dart() ).value() ) ) );
        };

        const topology::CombinatorialMapRestriction sides( bdry, keep_face_sides );
        const auto base_ptr = std::make_shared<const topology::CombinatorialMapRestriction>( bdry, keep_face_base, true );
        const auto& base = *base_ptr;
        const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

        if( eulerCharacteristic( base ) != 1 or eulerCharacteristic( target ) != 1 )
            throw std::runtime_error( "Input mesh is not the sweep of a disk-like topology." );

        const auto vertex_positions = [&sweep_input]( const topology::CombinatorialMap& map ) {
            const auto vertex_ids = indexingOrError( map, 0 );
            return [&sweep_input, vertex_ids]( const topology::Vertex& v ) -> Eigen::Vector3d {
                return sweep_input.mesh.points.at( vertex_ids( v ) );
            };
        };

        const topology::Edge start_edge = [&]() {
            std::optional<topology::Edge> out;
            iterateDartsWhile( sides, [&]( const topology::Dart& d ) {
                const auto phi2 = phi( bdry, 2, d );
                if( keep_face_base( phi2.value() ) )
                {
                    out.emplace( d );
                    return false;
                }
                return true;
            } );
            if( not out.has_value() ) throw std::runtime_error( "Unable to find tracing start edge" );
            return out.value();
        }();
        const reparam::Trace trace =
            reparam::traceBoundaryField( sides, start_edge, 0.5, sol, vertex_positions( sides ), false );

        LOG( log_level_set_based_tracing ) << "FINISHED TRACE\n\n";

        const std::vector<reparam::TraceLevelSetIntersection> intersections =
            reparam::levelSetIntersections( trace, sides, level_set_values );
        if( intersections.size() != level_set_values.size() )
            throw std::runtime_error( "Wrong number of intersections" );

        const auto vertex_ids = indexingOrError( map, 0 );

        const auto harmonic_func = [&]( const topology::Vertex& v ) { return sol( vertex_ids( v ) ); };

        const auto process_param =
            [&]( const std::shared_ptr<const topology::CombinatorialMap>& cmap,
                 const auto& positions,
                 const auto& thetas,
                 FoliationLeaf& leaf ) {
                const auto base_vert_ids = indexingOrError( *cmap, 0 );
                const std::map<size_t, double> thetas_by_id = [&]() {
                    std::map<size_t, double> out;
                    for( const auto& pr : thetas )
                    {
                        out.insert( { base_vert_ids( pr.first ), pr.second } );
                    }
                    return out;
                }();

                const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
                    if( boundaryAdjacent( *cmap, v ) )
                    {
                        const double theta = thetas_by_id.at( base_vert_ids( v ) );
                        return Eigen::Vector2d( cos( theta ), sin( theta ) );
                    }
                    else
                        return std::nullopt;
                };

                leaf.tutte = std::make_shared<const Eigen::MatrixX2d>(
                    reparam::tutteEmbedding( *cmap, positions, constraints_func, Laplace2dEdgeWeights::InverseLength ) );

                const auto atlas = std::make_shared<param::TriangleParametricAtlas>( cmap );
                const auto vert_positions = [tutte = *( leaf.tutte ), base_vert_ids]( const topology::Vertex& v ) {
                    return tutte.row( base_vert_ids( v ) );
                };
                leaf.tutte_mapping = std::make_shared<mapping::TriangleMeshCircleMapping>( atlas, vert_positions );
                leaf.space_mapping = std::make_shared<mapping::TriangleMeshMapping>( atlas, positions, 3 );
            };

        std::vector<FoliationLeaf> leaves;
        leaves.reserve( level_set_values.size() );

        { // Base level set
            const auto base_positions = vertex_positions( bdry );
            const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
                return map_face_ids( bdry.toUnderlyingCell( topology::Face( phi( bdry, 2, e.dart() ).value() ) ) );
            };
            const std::map<topology::Vertex, double> thetas =
                reparam::thetaValues( base, base_positions, face_ids_of_edge, intersections[0] );

            leaves.push_back( {} );
            process_param( base_ptr, base_positions, thetas, leaves.back() );
        }

        LOG( log_level_set_based_tracing ) << "FINISHED BASE\n\n";

        for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
        { // Midway level set
            leaves.push_back( {} );
            const auto level_set_cmap =
                std::make_shared<topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
            const auto v_pos = vertex_positions( map );
            const auto& level_set = *level_set_cmap;
            const auto level_set_positions = topology::levelSetVertexPositions( level_set, v_pos );
            const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
                return map_face_ids( level_set.underlyingCell( e ) );
            };
            const std::map<topology::Vertex, double> thetas =
                reparam::thetaValues( level_set, level_set_positions, face_ids_of_edge, intersections[level_ii] );

            const auto level_set_tri =
                std::make_shared<topology::DelaunayTriangulation>( level_set_cmap, level_set_positions );
            const auto tri_positions =
                topology::delaunayTriangulationVertexPositions( *level_set_tri, level_set_positions );

            process_param( level_set_tri, tri_positions, thetas, leaves.back() );

            LOG( log_level_set_based_tracing ) << "FINISHED LEVEL " << ( level_ii + 1 ) << std::endl << std::endl;
        }

        { // target level set
            leaves.push_back( {} );
            const auto target_positions = vertex_positions( bdry );
            const auto reversed_cmap = std::make_shared<topology::ReversedCombinatorialMap>( target );
            const auto& rev_map = *reversed_cmap;
            const auto rev_positions = reversedVertexPositions( rev_map, target_positions );
            const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
                return map_face_ids( bdry.toUnderlyingCell(
                    topology::Face( phi( bdry, 2, rev_map.toUnderlyingCell( e ).dart() ).value() ) ) );
            };
            const std::map<topology::Vertex, double> thetas =
                reparam::thetaValues( rev_map, rev_positions, face_ids_of_edge, intersections.back() );

            process_param( reversed_cmap, rev_positions, thetas, leaves.back() );
        }
        LOG( log_level_set_based_tracing ) << "FINISHED TARGET\n\n";

        callback( leaves );
    }

    Eigen::MatrixXd fitLinearMeshToLeaves(
        const basis::SplineSpace& ss,
        const std::vector<reparam::FoliationLeaf>& leaves,
        const std::function<std::pair<Eigen::Vector2d, size_t>( const topology::Vertex& )>& tutte_points,
        const std::optional<std::function<Eigen::Vector3d( const topology::Vertex& )>>& first_level_points )
    {
        const auto sample_at = [&]( const topology::Vertex& v ) -> Eigen::Vector3d {
            // Get the 2d vertex and the level set that correspond to this vertex.
            const auto [tutte_pt, leaf_ii] = tutte_points( v );
            if( leaf_ii == 0 and first_level_points ) return first_level_points.value()( v );

            const Eigen::Vector3d field_pt = [&]() -> Eigen::Vector3d {
                const auto param_pt = leaves.at( leaf_ii ).tutte_mapping->maybeInverse( tutte_pt );

                if( param_pt.has_value() )
                    return leaves.at( leaf_ii ).space_mapping->evaluate( param_pt.value().first,
                                                                         param_pt.value().second );
                else
                    throw std::runtime_error( "Failed to find a tutte point inverse!" );
            }();

            return field_pt;
        };

        Eigen::MatrixXd fit_cpts = Eigen::MatrixXd::Zero( cellCount( ss.basisComplex().parametricAtlas().cmap(), 0 ), 3 );

        iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 0, [&]( const topology::Vertex& v ) {
            const auto conn = ss.connectivity( v );
            if( conn.size() != 1 )
            {
                return true;
            }

            fit_cpts.row( conn.front() ) = sample_at( v );
            return true;
        } );

        return fit_cpts;
    }
} // namespace reparam