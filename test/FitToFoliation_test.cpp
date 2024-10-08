#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
#include <TriangleParametricAtlas.hpp>
#include <TriangleMeshCircleMapping.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <CombinatorialMapRestriction.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <SimplexUtilities.hpp>
#include <Laplace.hpp>
#include <Tracing.hpp>
#include <Foliation.hpp>
#include <LevelSetCMap.hpp>
#include <ReversedCombinatorialMap.hpp>
#include <DelaunayTriangulation.hpp>
#include <VTKOutput.hpp>
#include <AbaqusInput.hpp>
#include <iomanip>
#include <memory>
#include <SplineSpaceEvaluator.hpp>
#include <LeastSquaresFitting.hpp>
#include <TPSplineSpace.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ranges>

struct FoliationLeaf
{
    // For the middle leaves
    std::shared_ptr<topology::LevelSetCMap> level_set_cmap;
    std::shared_ptr<topology::DelaunayTriangulation> level_set_tri;

    // For the target surface
    std::shared_ptr<topology::ReversedCombinatorialMap> reversed_cmap;

    // For all leaves
    std::shared_ptr<const Eigen::MatrixX2d> tutte;
    std::shared_ptr<param::TriangleParametricAtlas> atlas;
    std::shared_ptr<mapping::TriangleMeshCircleMapping> circle_mapping;
    std::shared_ptr<mapping::TriangleMeshMapping> space_mapping;
};

void levelSetBasedTracing( const SweepInput& sweep_input,
                           const std::vector<double> level_set_values,
                           const bool log,
                           const std::optional<std::string> output_prefix,
                           const std::function<void(const std::vector<FoliationLeaf>&)>& callback )
{
    const topology::TetMeshCombinatorialMap map( sweep_input.mesh );
    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd sol = reparam::sweepEmbedding( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    if( log ) std::cout << "FINISHED LAPLACE\n\n";

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
    const topology::CombinatorialMapRestriction base( bdry, keep_face_base, true );
    const topology::CombinatorialMapRestriction target( bdry, keep_face_target, true );

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

    if( log ) std::cout << "FINISHED TRACE\n\n";

    const std::vector<reparam::TraceLevelSetIntersection> intersections =
        reparam::levelSetIntersections( trace, sides, level_set_values );
    REQUIRE( intersections.size() == level_set_values.size() );

    const auto vertex_ids = indexingOrError( map, 0 );

    const auto harmonic_func = [&]( const topology::Vertex& v ) { return sol( vertex_ids( v ) ); };

    SimplicialComplex level_sets;

    const auto process_param =
        [&]( const topology::CombinatorialMap& cmap, const auto& positions, const auto& thetas, FoliationLeaf& leaf ) {
            const auto base_vert_ids = indexingOrError( cmap, 0 );
            const std::map<size_t, double> thetas_by_id = [&]() {
                std::map<size_t, double> out;
                for( const auto& pr : thetas )
                {
                    out.insert( { base_vert_ids( pr.first ), pr.second } );
                }
                return out;
            }();

            const auto constraints_func = [&]( const topology::Vertex& v ) -> std::optional<Eigen::Vector2d> {
                if( boundaryAdjacent( cmap, v ) )
                {
                    const double theta = thetas_by_id.at( base_vert_ids( v ) );
                    return Eigen::Vector2d( cos( theta ), sin( theta ) );
                }
                else
                    return std::nullopt;
            };

            leaf.tutte = std::make_unique<const Eigen::MatrixX2d>(
                reparam::tutteEmbedding( cmap, positions, constraints_func, true ) );

            leaf.atlas = std::make_unique<param::TriangleParametricAtlas>( cmap );
            const auto vert_positions = [tutte = *( leaf.tutte ), base_vert_ids]( const topology::Vertex& v ) {
                return tutte.row( base_vert_ids( v ) );
            };
            leaf.circle_mapping = std::make_unique<mapping::TriangleMeshCircleMapping>( *leaf.atlas, vert_positions );
            leaf.space_mapping = std::make_unique<mapping::TriangleMeshMapping>( *leaf.atlas, positions, 3 );

            if( output_prefix )
                iterateCellsWhile( cmap, 2, [&]( const topology::Face& f ) {
                    addTriangleNoDuplicateChecking( level_sets, triangleOfFace<3>( cmap, positions, f ) );
                    return true;
                } );
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
        process_param( base, base_positions, thetas, leaves.back() );
    }

    if( log ) std::cout << "FINISHED BASE\n\n";

    for( size_t level_ii = 1; level_ii < level_set_values.size() - 1; level_ii++ )
    { // Midway level set
        leaves.push_back( {} );
        leaves.back().level_set_cmap =
            std::make_unique<topology::LevelSetCMap>( map, harmonic_func, level_set_values[level_ii] );
        const auto v_pos = vertex_positions( map );
        const auto& level_set = *leaves.back().level_set_cmap;
        const auto level_set_positions = topology::levelSetVertexPositions( level_set, v_pos );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
            return map_face_ids( level_set.underlyingCell( e ) );
        };
        const std::map<topology::Vertex, double> thetas =
            reparam::thetaValues( level_set, level_set_positions, face_ids_of_edge, intersections[level_ii] );

        leaves.back().level_set_tri =
            std::make_unique<topology::DelaunayTriangulation>( level_set, level_set_positions );
        const auto tri_positions =
            topology::delaunayTriangulationVertexPositions( *leaves.back().level_set_tri, level_set_positions );

        process_param( *leaves.back().level_set_tri, tri_positions, thetas, leaves.back() );

        if( log ) std::cout << "FINISHED LEVEL " << ( level_ii + 1 ) << std::endl << std::endl;
    }

    { // target level set
        leaves.push_back( {} );
        const auto target_positions = vertex_positions( bdry );
        leaves.back().reversed_cmap = std::make_unique<topology::ReversedCombinatorialMap>( target );
        const auto& rev_map = *leaves.back().reversed_cmap;
        const auto rev_positions = reversedVertexPositions( rev_map, target_positions );
        const auto face_ids_of_edge = [&]( const topology::Edge& e ) {
            return map_face_ids( bdry.toUnderlyingCell(
                topology::Face( phi( bdry, 2, rev_map.toUnderlyingCell( e ).dart() ).value() ) ) );
        };
        const std::map<topology::Vertex, double> thetas =
            reparam::thetaValues( rev_map, rev_positions, face_ids_of_edge, intersections.back() );

        process_param( rev_map, rev_positions, thetas, leaves.back() );
    }
    if( log ) std::cout << "FINISHED TARGET\n\n";

    if( output_prefix )
    {
        // Draw some lines!
        SimplicialComplex param_out;
        param_out.points.reserve( 8 * 6 * level_set_values.size() );
        param_out.simplices.reserve( 8 * 6 * level_set_values.size() );
        for( const double radius : { 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 } )
        {
            for( const double theta : { 0, 45, 90, 135, 180, 225, 270, 315 } )
            {
                const Eigen::Vector2d circle_pt( radius * cos( theta * std::numbers::pi / 180 ),
                                                radius * sin( theta * std::numbers::pi / 180 ) );
                const size_t offset = param_out.points.size();
                size_t i = 0;
                for( const auto& leaf : leaves )
                {
                    const auto& param_pt = leaf.circle_mapping->maybeInverse( circle_pt );
                    CHECK( param_pt.has_value() );
                    if( not param_pt.has_value() )
                    {
                        std::cerr << "NO VALUE r: " << radius << " t: " << theta << " i: " << i << std::endl;
                        break;
                    }
                    if( output_prefix )
                    {
                        const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                        param_out.points.push_back( space_pt );
                    }
                    i++;
                }

                if( output_prefix )
                {
                    if( i == 0 ) continue;
                    for( size_t simplex_ii = 0; simplex_ii < i - 1; simplex_ii++ )
                    {
                        param_out.simplices.push_back( Simplex( offset + simplex_ii, offset + simplex_ii + 1 ) );
                    }
                }
            }
        }

        io::VTKOutputObject output( level_sets );
        io::outputSimplicialFieldToVTK( output, output_prefix.value() + "_level_sets.vtu" );

        io::VTKOutputObject output2( param_out );
        io::outputSimplicialFieldToVTK( output2, output_prefix.value() + "_foliation_param.vtu" );
    }

    callback( leaves );
}

void levelSetBasedTracing( const SweepInput& sweep_input,
                           const size_t n_levels,
                           const bool log,
                           const std::optional<std::string> output_prefix,
                           const std::function<void( const std::vector<FoliationLeaf>& )>& callback )
{
    const std::vector<double> level_set_values = [&]( const size_t n_levels ) {
        std::vector<double> out;
        out.reserve( n_levels );
        const double diff = 1.0 / ( n_levels - 1 );
        for( size_t i = 0; i < n_levels; i++ )
        {
            out.push_back( i * diff );
        }
        return out;
    }( n_levels );

    levelSetBasedTracing( sweep_input, level_set_values, log, output_prefix, callback );
}


std::vector<double> linspace( const double left_val, const double right_val, const size_t n_levels ) {
    std::vector<double> out;
    out.reserve( n_levels );
    const double diff = ( right_val - left_val ) / ( n_levels - 1 );
    for( size_t i = 0; i < n_levels; i++ )
    {
        out.push_back( left_val + i * diff );
    }
    return out;
}

std::vector<double> concatenate( const std::vector<double>& first, const std::vector<double>& second ) {
    std::vector<double> out;
    out.reserve( first.size() + second.size() );

    out.insert( out.end(), first.begin(), first.end() );
    out.insert( out.end(), second.begin(), second.end() );

    return out;
}

std::vector<std::pair<topology::Cell, param::ParentPoint>> parentPointsOfParamPoints(
    const std::vector<double>& level_set_values, const size_t n_elems_u, const double param_tol )
{
    std::vector<std::pair<topology::Cell, param::ParentPoint>> out;
    out.reserve( level_set_values.size() );

    const double elem_extent = 1.0 / n_elems_u;

    const param::ParentDomain pd = param::cubeDomain( 1 );

    for( const double value : level_set_values )
    {
        if( value < -param_tol or value > 1 + param_tol )
            throw std::out_of_range( "Value outside of [0, 1] range" );

        const double clamped_value = std::max( 0.0, std::min( 1.0, value ) );

        const size_t elem = [&](){
            const size_t elem = static_cast<size_t>( std::floor( clamped_value / elem_extent ) );
            return elem == n_elems_u ? elem - 1 : elem;
        }();

        const double elem_start = elem * elem_extent;
        const double u = ( clamped_value - elem_start ) / elem_extent;

        out.push_back( { topology::Edge( elem ),
                         param::ParentPoint( pd,
                                             Vector1d( u ),
                                             { util::equals( u, 1.0, param_tol ),
                                               util::equals( u, 0.0, param_tol ) } ) } );
    }

    return out;
}

std::vector<std::pair<topology::Cell, param::ParentPoint>> parentPointsOfParamPoints(
    const std::vector<double>& values, const param::ParametricAtlas1d& pa, const double param_tol )
{
    std::vector<std::pair<topology::Cell, param::ParentPoint>> out;
    out.reserve( values.size() );

    size_t value_idx = 0;
    double current_position = values.front();
    const double factor = ( values.back() - values.front() ) / pa.totalLength();

    const param::ParentDomain pd = param::cubeDomain( 1 );

    iterateCellsWhile( pa.cmap(), 1, [&]( const topology::Edge& c ) {
        const double interval_length = pa.parametricLengths( c )( 0 ) * factor;
        const double next_position = current_position + interval_length;

        // Process all values that fall within the current interval
        for( ; value_idx < values.size() and values.at( value_idx ) <= next_position; value_idx++ )
        {
            const auto [relative_pos, zerovec] = [&]() -> std::pair<double, param::BaryCoordIsZeroVec> {
                if( values.at( value_idx ) <= current_position + param_tol )
                {
                    return { 0.0, { false, true } };
                }
                if( values.at( value_idx ) >= next_position - param_tol )
                {
                    return { 1.0, { true, false } };
                }
                return { ( values.at( value_idx ) - current_position ) / interval_length, { false, false } };
            }();

            out.emplace_back( c, param::ParentPoint( pd, Vector1d( relative_pos ), zerovec ) );
        }

        current_position = next_position;
        return true;
    } );

    for( ; value_idx < values.size(); value_idx++ )
    {
        if( values.at( value_idx ) <= current_position + param_tol )
            out.emplace_back( topology::Edge( pa.cmap().maxDartId() ), param::ParentPoint( pd, Vector1d( 1.0 ), { true, false } ) );
        else
            throw std::runtime_error( "Level set values outside of parametric domain" );
    }

    std::cout << out << std::endl;

    return out;
}

Eigen::Vector2d toUnitSquare( const topology::TPCombinatorialMap& cmap, const topology::Face& f, const Eigen::Vector2d& pt )
{
    // FIXME: Assumes unit parametric lengths
    const auto [d1, d2, _] = cmap.unflatten( f.dart() );
    return Eigen::Vector2d( ( (double)d1.id() + pt( 0 ) ) / (double)cellCount( cmap.sourceCMap(), 1 ), ( (double)d2.id() + pt( 1 ) ) / (double)cellCount( cmap.lineCMap(), 1 ) );
}

// transforms [0,1]x[0,1] to the unit circle.  Found with a public domain license at
// https://github.com/erich666/jgt-code/blob/master/Volume_02/Number_3/Shirley1997/disk.cc
Eigen::Vector2d toUnitDisk1( const Eigen::Vector2d& square_coords )
{
    double phi, r;
    const double a = 2 * square_coords( 0 ) - 1;
    const double b = 2 * square_coords( 1 ) - 1;
    if( a == 0 and b == 0 )
    {
        r = 0;
        phi = 0;
    }
    if( a * a > b * b )
    {
        r = a;
        phi = ( std::numbers::pi / 4 ) * ( b / a );
    }
    else
    {
        r = b;
        phi = ( std::numbers::pi / 2 ) - ( std::numbers::pi / 4 ) * ( a / b );
    }
    return Eigen::Vector2d( r * cos( phi ), r * sin( phi ) );
}

Eigen::Vector2d toUnitDisk( const Eigen::Vector2d& square_coords )
{
    const double a = 2 * square_coords( 0 ) - 1;
    const double b = 2 * square_coords( 1 ) - 1;
    if( a == 0 and b == 0 ) return Eigen::Vector2d( a, b );
    const double asq = a*a;
    const double bsq = b*b;
    const double scaling = sqrt( asq + bsq - asq * bsq ) / sqrt( asq + bsq );
    return scaling * Eigen::Vector2d( a, b );
}

Eigen::Vector2d multiPatchToUnitDisk( const size_t patch_ii, const Eigen::Vector2d& square_coords )
{
    constexpr double a = std::numbers::sqrt2 / 4;
    using std::numbers::pi;
    const double& s = square_coords( 0 );
    const double& t = square_coords( 1 );
    switch( patch_ii )
    {
        case 0:
            return a * 2 * square_coords - Eigen::Vector2d::Constant( a );
        case 1:
            return Eigen::Vector2d( ( 1 - s ) * a + s * cos( 0.25 * ( 2 * t - 1 ) * pi ),
                                    ( 2 * t - 1 ) * a * ( 1 - s ) + s * sin( 0.25 * ( 2 * t - 1 ) * pi ) );
        case 2:
            return Eigen::Vector2d( ( 1 - 2 * t ) * a * ( 1 - s ) + s * cos( 0.25 * pi + t * 0.5 * pi ),
                                    ( 1 - s ) * a + s * sin( 0.25 * pi + t * 0.5 * pi ) );
        case 3:
            return Eigen::Vector2d( ( s - 1 ) * a + s * cos( 0.75 * pi + t * 0.5 * pi ),
                                    ( 1 - 2 * t ) * a * ( 1 - s ) + s * sin( 0.75 * pi + t * 0.5 * pi ) );
        case 4:
            return Eigen::Vector2d( ( 2 * t - 1 ) * a * ( 1 - s ) + s * cos( 1.25 * pi + t * 0.5 * pi ),
                                    ( s - 1 ) * a + s * sin( 1.25 * pi + t * 0.5 * pi ) );
        default:
            throw std::runtime_error( "Bad patch id" );
    }
}

void fitToPringlesSinglePatch( const SweepInput& sweep_input,
                               const std::vector<double>& level_set_values,
                               const bool log_progress,
                               const std::string& output_prefix,
                               const size_t degree,
                               const size_t n_elems_u,
                               const size_t n_elems_st )
{
    std::vector<double> level_set_values2;
    level_set_values2.push_back( 0 );
    level_set_values2.insert( level_set_values2.end(), level_set_values.begin(), level_set_values.end() );
    level_set_values2.push_back( 1.0 );
    levelSetBasedTracing( sweep_input, level_set_values2, log_progress, std::nullopt, [&]( const std::vector<FoliationLeaf>& leaves2 ) {
        std::vector<FoliationLeaf> leaves( std::next( leaves2.begin() ), leaves2.end() );
        leaves.pop_back();
        const std::vector<std::pair<topology::Cell, param::ParentPoint>> ppt_u = parentPointsOfParamPoints( level_set_values, n_elems_u, 1e-9 );

        const basis::KnotVector kv1 = basis::integerKnotsWithNElems( n_elems_st, degree );
        const basis::KnotVector kv2 = basis::integerKnotsWithNElems( n_elems_u, degree );

        // TODO: use multipatch
        const basis::TPSplineSpace vol_ss = basis::buildBSpline( {kv1, kv1, kv2}, {degree, degree, degree} );
        const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( vol_ss.source() );
        const topology::TPCombinatorialMap& vol_cmap = vol_ss.basisComplex().parametricAtlas().cmap();

        const size_t n_points = cellCount( source_ss.basisComplex().parametricAtlas().cmap(), 2 ) * 4 * leaves.size();

        const SmallVector<double, 2> source_points{ 0.1, 0.9 };

        const param::ParentDomain pd_3d = param::cubeDomain( 3 );

        std::cout << "About to fit\n";

        SimplicialComplex fitting_points;

        eval::SplineSpaceEvaluator evaler( vol_ss, 0 );
        const auto fit_cpts = fitting::leastSquaresFitting(
            evaler,
            n_points,
            3,
            [&]( const std::function<void( const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd& )>&
                    callback ) {
                        // std::cout << "Inside callback!\n";
                for( size_t leaf_ii = 0; leaf_ii < leaves.size(); leaf_ii++ )
                {
                    // std::cout << "Leaf " << leaf_ii << std::endl;
                    iterateCellsWhile(
                        source_ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const topology::Cell& f ) {
                            // std::cout << "Source cell " << f << std::endl;
                            util::iterateTensorProduct(
                                { source_points.size(), source_points.size() }, [&]( const util::IndexVec& indices ) {
                                    const Eigen::Vector3d pt( source_points.at( indices.at( 0 ) ),
                                                            source_points.at( indices.at( 1 ) ),
                                                            ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
                                    // std::cout << "Point: " << pt.transpose() << std::endl;
                                    const param::ParentPoint vol_ppt(
                                        pd_3d,
                                        pt,
                                        { false,
                                        false,
                                        false,
                                        false,
                                        ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 0 ),
                                        ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 1 ) } );

                                    const topology::Volume cell(
                                        vol_cmap.flatten( f.dart(),
                                                        ppt_u.at( leaf_ii ).first.dart(),
                                                        topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );

                                    const Eigen::Vector2d circle_pt = toUnitDisk( toUnitSquare( source_ss.basisComplex().parametricAtlas().cmap(), f, pt.head( 2 ) ) );

                                    // std::cout << "About to invert circle point " << circle_pt.transpose() << std::endl;
                                    // std::cout << "Circle mapping spatial dims: " << leaves.at( leaf_ii ).circle_mapping->spatialDim() << std::endl;
                                    const auto param_pt = leaves.at( leaf_ii ).circle_mapping->maybeInverse( circle_pt );
                                    // std::cout << "Inverted it!" << std::endl;
                                    CHECK( param_pt.has_value() );
                                    if( not param_pt.has_value() )
                                    {
                                        std::cerr << "NO VALUE" << std::endl;
                                        pauseDebugger();
                                    }
                                    const Eigen::Vector3d field_pt = leaves.at( leaf_ii ).space_mapping->evaluate( param_pt.value().first, param_pt.value().second );

                                    fitting_points.simplices.emplace_back( fitting_points.points.size() );
                                    fitting_points.points.push_back( field_pt );

                                    // std::cout << "Call back on " << cell << " " << vol_ppt << " " << field_pt.transpose() << std::endl;
                                    callback( cell, vol_ppt, field_pt );
                                } );
                            return true;
                        } );
                }
            } );

        io::VTKOutputObject fitting_points_output( fitting_points );
        io::outputSimplicialFieldToVTK( fitting_points_output, output_prefix + "fitting_points.vtu" );

        std::cout << "Finished fit\n";

        std::cout << "Vol cells: " << cellCount( vol_cmap, 3 ) << std::endl;

        io::outputBezierMeshToVTK( vol_ss,
                                    fit_cpts,
                                    "fit_to_" + output_prefix + ".vtu" );
    } );
}

void fitToPringles5Patch( const SweepInput& sweep_input,
                          const std::vector<double>& level_set_values,
                          const bool log_progress,
                          const std::string& output_prefix,
                          const size_t degree,
                          const basis::KnotVector& kv_u,
                          const size_t n_elems_st )
{
    using namespace topology;
    std::vector<double> level_set_values2;
    const bool add_front = not util::equals( level_set_values.front(), 0.0, 1e-9 );
    const bool add_back = not util::equals( level_set_values.back(), 1.0, 1e-9 );
    if( add_front ) level_set_values2.push_back( 0 );
    level_set_values2.insert( level_set_values2.end(), level_set_values.begin(), level_set_values.end() );
    if( add_back ) level_set_values2.push_back( 1.0 );
    levelSetBasedTracing( sweep_input, level_set_values2, log_progress, std::nullopt, [&]( const std::vector<FoliationLeaf>& leaves2 ) {
        const std::vector<FoliationLeaf> leaves( add_front ? std::next( leaves2.begin() ) : leaves2.begin(), add_back ? std::prev( leaves2.end() ) : leaves2.end() );
        const basis::KnotVector kv_st = basis::integerKnotsWithNElems( n_elems_st, degree );

        const std::shared_ptr<const basis::TPSplineSpace> TP_ss = std::make_shared<const basis::TPSplineSpace>( basis::buildBSpline( {kv_st, kv_st, kv_u}, {degree, degree, degree} ) );
        const auto& vol_cmap = TP_ss->basisComplex().parametricAtlas().cmap();
        const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( TP_ss->source() );
        const std::vector<std::pair<topology::Cell, param::ParentPoint>> ppt_u = parentPointsOfParamPoints( level_set_values, TP_ss->line().basisComplex().parametricAtlas(), 1e-9 );

        const basis::MultiPatchSplineSpace mp_ss = basis::buildMultiPatchSplineSpace(
            std::vector<std::shared_ptr<const basis::TPSplineSpace>>( 5, TP_ss ),
            { // FIXME: only valid for n_elems_st = 2
                { { 0, Dart( 31 ) }, { 1, Dart( 19 ) } },
                { { 0, Dart( 85 ) }, { 2, Dart( 19 ) } },
                { { 0, Dart( 67 ) }, { 3, Dart( 19 ) } },
                { { 0, Dart( 1 ) }, { 4, Dart( 19 ) } },
                { { 1, Dart( 61 ) }, { 2, Dart( 1 ) } },
                { { 2, Dart( 61 ) }, { 3, Dart( 1 ) } },
                { { 3, Dart( 61 ) }, { 4, Dart( 1 ) } },
                { { 4, Dart( 61 ) }, { 1, Dart( 1 ) } },
            } );

        const size_t n_points = cellCount( source_ss.basisComplex().parametricAtlas().cmap(), 2 ) * 4 * 5 * leaves.size();

        const SmallVector<double, 2> source_points{ 0.1, 0.9 };

        const param::ParentDomain pd_3d = param::cubeDomain( 3 );

        std::cout << "About to fit\n";

        SimplicialComplex fitting_points;

        eval::SplineSpaceEvaluator evaler( mp_ss, 0 );
        const auto fit_cpts = fitting::leastSquaresFitting(
            evaler,
            n_points,
            3,
            [&]( const std::function<void( const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd& )>&
                    callback ) {
                        // std::cout << "Inside callback!\n";
                for( size_t leaf_ii = 0; leaf_ii < leaves.size(); leaf_ii++ )
                {
                    for( size_t patch_ii = 0; patch_ii < 5; patch_ii++ )
                    {
                        // std::cout << "Leaf " << leaf_ii << std::endl;
                        iterateCellsWhile(
                            source_ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const topology::Cell& f ) {
                                // std::cout << "Source cell " << f << std::endl;
                                util::iterateTensorProduct(
                                    { source_points.size(), source_points.size() }, [&]( const util::IndexVec& indices ) {
                                        const Eigen::Vector3d pt( source_points.at( indices.at( 0 ) ),
                                                                source_points.at( indices.at( 1 ) ),
                                                                ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
                                        // std::cout << "Point: " << pt.transpose() << std::endl;
                                        const param::ParentPoint vol_ppt(
                                            pd_3d,
                                            pt,
                                            { false,
                                            false,
                                            false,
                                            false,
                                            ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 0 ),
                                            ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 1 ) } );

                                        const topology::Volume patch_cell(
                                            vol_cmap.flatten( f.dart(),
                                                            ppt_u.at( leaf_ii ).first.dart(),
                                                            topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );

                                        const topology::Volume cell( mp_ss.basisComplex().parametricAtlas().cmap().toGlobalDart( patch_ii, patch_cell.dart() ) );

                                        const Eigen::Vector2d circle_pt = multiPatchToUnitDisk(
                                            patch_ii,
                                            toUnitSquare(
                                                source_ss.basisComplex().parametricAtlas().cmap(), f, pt.head( 2 ) ) );

                                        // std::cout << "About to invert circle point " << circle_pt.transpose() << std::endl;
                                        // std::cout << "Circle mapping spatial dims: " << leaves.at( leaf_ii ).circle_mapping->spatialDim() << std::endl;
                                        const auto param_pt = leaves.at( leaf_ii ).circle_mapping->maybeInverse( circle_pt );
                                        // std::cout << "Inverted it!" << std::endl;
                                        CHECK( param_pt.has_value() );
                                        if( not param_pt.has_value() )
                                        {
                                            std::cerr << "NO VALUE" << std::endl;
                                            pauseDebugger();
                                        }
                                        const Eigen::Vector3d field_pt = leaves.at( leaf_ii ).space_mapping->evaluate( param_pt.value().first, param_pt.value().second );

                                        fitting_points.simplices.emplace_back( fitting_points.points.size() );
                                        fitting_points.points.push_back( field_pt );

                                        // std::cout << "Call back on " << cell << " " << vol_ppt << " " << field_pt.transpose() << std::endl;
                                        callback( cell, vol_ppt, field_pt );
                                    } );
                                return true;
                            } );
                    }
                }
            } );

        io::VTKOutputObject fitting_points_output( fitting_points );
        io::outputSimplicialFieldToVTK( fitting_points_output, output_prefix + "fitting_points.vtu" );

        io::outputBezierMeshToVTK( mp_ss,
                                    fit_cpts,
                                    "fit_to_" + output_prefix + "_multi_patch.vtu" );
    } );
}

void fitToPringles5Patch( const SweepInput& sweep_input,
                          const std::vector<double>& level_set_values,
                          const bool log_progress,
                          const std::string& output_prefix,
                          const size_t degree,
                          const size_t& n_elems_u,
                          const size_t n_elems_st )
{
    const basis::KnotVector kv_u = basis::integerKnotsWithNElems( n_elems_u, degree );
    return fitToPringles5Patch( sweep_input, level_set_values, log_progress, output_prefix, degree, kv_u, n_elems_st );
}


TEST_CASE( "Level set parameterization of left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );

    // level set every 0.02 until 0.75, then ~30 levels between 0.75 and 0.85, then every 0.2 until 1.0.
    // const std::vector<double> level_set_values =
    //     concatenate(
    //         concatenate( linspace( 0, 0.78, 40 ), linspace( 0.8, 0.819, 20 ) ),
    //         concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 31 ) ), linspace( 0.86, 1.0, 13 ) ) );
    // const std::vector<double> level_set_values =
    //     concatenate(
    //         concatenate( linspace( 0, 0.78, 20 ), linspace( 0.8, 0.819, 10 ) ),
    //         concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 15 ) ), linspace( 0.86, 1.0, 7 ) ) );
    const std::vector<double> level_set_values = linspace( 0, 1.0, 80 );
    const bool log_progress = false;
    const std::string output_prefix = "ventricle";

    // fitToPringlesSinglePatch( sweep_input, level_set_values, log_progress, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, log_progress, output_prefix, 2, 20, 2 );
}

TEST_CASE( "Level set parameterization of part of left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );

    // level set every 0.02 until 0.75, then ~30 levels between 0.75 and 0.85, then every 0.2 until 1.0.
    // const std::vector<double> level_set_values =
    //     concatenate(
    //         concatenate( linspace( 0, 0.78, 40 ), linspace( 0.8, 0.819, 20 ) ),
    //         concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 31 ) ), linspace( 0.86, 1.0, 13 ) ) );
    // const std::vector<double> level_set_values =
    //     concatenate(
    //         concatenate( linspace( 0, 0.78, 20 ), linspace( 0.8, 0.819, 10 ) ),
    //         concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 15 ) ), linspace( 0.86, 1.0, 7 ) ) );
    // const std::vector<double> level_set_values = linspace( 0, 1.0, 80 );
    // const std::vector<double> level_set_values = linspace( 0.75, 0.85, 80 );
    const std::vector<double> level_set_values =
        concatenate( concatenate( linspace( 0.75, 0.8, 4 ), linspace( 0.8, 0.818, 14 ) ),
                     concatenate( linspace( 0.81825, 0.82275, 19 ), linspace( 0.823, 0.85, 31 ) ) );
    const bool log_progress = false;
    const std::string output_prefix = "part_ventricle";

    // fitToPringlesSinglePatch( sweep_input, level_set_values, log_progress, output_prefix, 2, 20, 4 );

    // const basis::KnotVector kv_u( {
    //     0, 0, 0, 6, 10, 12, 12.5, 13, 13.2, 13.4, 13.6, 13.8, 14, 14.2, 14.4, 14.6, 14.8, 15, 15.5, 16, 17, 18, 20, 20, 20 
    // }, 1e-9 );
    const basis::KnotVector kv_u( {
        0, 0, 0, 6, 12, 12.7, 13, 13.25, 13.5, 13.75, 14, 14.25, 14.5, 14.75, 15, 15.4, 16, 17, 18, 19, 20, 20, 20 
    }, 1e-9 );
    fitToPringles5Patch( sweep_input, level_set_values, log_progress, output_prefix, 2, kv_u, 2 );
}

TEST_CASE( "Level set parameterization of hook" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/hook.inp", "Surface12", "Surface10" );

    const bool log_progress = false;
    const std::string output_prefix = "hook";
    const std::vector<double> level_set_values = linspace( 0, 1.0, 30 );

    // fitToPringlesSinglePatch( sweep_input, level_set_values, log_progress, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, log_progress, output_prefix, 2, 20, 2 );
}

// TEST_CASE( "Level set parameterization of cube" )
// {
//     const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();

//     const bool log_progress = false;
//     levelSetBasedTracing( sweep_input, 3, log_progress );
// }