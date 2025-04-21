#include <SweepAPIMethods.hpp>
#include <SweepInput.hpp>
#include <Foliation.hpp>
#include <VTKOutput.hpp>
#include <TriangleMeshMapping.hpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMapMethods.hpp>

#include <KnotVector.hpp>
#include <CommonUtils.hpp>
#include <TPSplineSpace.hpp>
#include <MultiPatchSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <iostream>

namespace api
{
    void outputLevelSetsAndTraces( const Sweep& sweep,
                                   const std::vector<double>& level_set_values,
                                   const std::vector<Eigen::Vector2d>& trace_points,
                                   const std::string& output_prefix )
    {
        if( sweep.source.size() == 0 or sweep.target.size() == 0 )
            throw std::invalid_argument( "No source or target surface specified" );
        if( sweep.mesh.points.size() == 0 or sweep.mesh.simplices.size() == 0 )
            throw std::invalid_argument( "No tet mesh supplied" );
        if( level_set_values.size() < 2 )
            throw std::invalid_argument( "Insufficient level sets specified; please include at least two values." );

        const SweepInput sweep_input = [&sweep]() {
            const auto& m = sweep.mesh;
            std::vector<bool> zero_bcs( m.points.size(), false );
            std::vector<bool> one_bcs( m.points.size(), false );
            for( const VertexId::Type& vid : sweep.source ) zero_bcs.at( vid ) = true;
            for( const VertexId::Type& vid : sweep.target ) one_bcs.at( vid ) = true;
            return SweepInput( m, zero_bcs, one_bcs );
        }();

        reparam::levelSetBasedTracing( sweep_input, level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
            SimplicialComplex level_sets;
            SimplicialComplex param_out;

            for( const auto& leaf : leaves )
            {
                iterateCellsWhile( leaf.space_mapping->parametricAtlas().cmap(), 2, [&]( const topology::Face& f ) {
                    addTriangleNoDuplicateChecking( level_sets,
                                                    triangleOfFace<3>( leaf.space_mapping->parametricAtlas().cmap(),
                                                                    leaf.space_mapping->vertPositions(),
                                                                    f ) );
                    return true;
                } );
            }

            // Draw some lines!
            param_out.points.reserve( trace_points.size() * level_set_values.size() );
            param_out.simplices.reserve( trace_points.size() * level_set_values.size() );

            for( const auto& circle_pt : trace_points )
            {
                const size_t offset = param_out.points.size();
                size_t i = 0;
                for( const auto& leaf : leaves )
                {
                    const auto& param_pt = leaf.tutte_mapping->maybeInverse( circle_pt );
                    if( not param_pt.has_value() )
                    {
                        std::cerr << "NO VALUE pt: " << circle_pt.transpose() << " level: " << i << std::endl;
                        break;
                    }
                    const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                    param_out.points.push_back( space_pt );
                    i++;
                }

                if( i == 0 ) continue;
                for( size_t simplex_ii = 0; simplex_ii < i - 1; simplex_ii++ )
                {
                    param_out.simplices.push_back( Simplex( offset + simplex_ii, offset + simplex_ii + 1 ) );
                }
            }

            io::VTKOutputObject output( level_sets );
            io::outputSimplicialFieldToVTK( output, output_prefix + "_level_sets.vtu" );

            io::VTKOutputObject output2( param_out );
            io::outputSimplicialFieldToVTK( output2, output_prefix + "_traces.vtu" );

            io::VTKOutputObject output3( sweep.mesh );
            io::outputSimplicialFieldToVTK( output3, output_prefix + "_tet_mesh.vtu" );
        } );
    }

    HexMesh fitSinglePatchHexMeshToSweep( const api::Sweep& sweep, const size_t n_elems_st, const std::vector<double>& u_values, const bool debug )
    {
        if( sweep.source.size() == 0 or sweep.target.size() == 0 )
            throw std::invalid_argument( "No source or target surface specified" );
        if( sweep.mesh.points.size() == 0 or sweep.mesh.simplices.size() == 0 )
            throw std::invalid_argument( "No tet mesh supplied" );
        if( u_values.size() < 2 )
            throw std::invalid_argument( "Insufficient u coordinates specified; please include at least two values." );
        if( n_elems_st < 1 )
            throw std::invalid_argument( "Invalid number of elements in the s and t directions specified" );
        if( u_values.front() != 0.0 or u_values.back() != 1.0 )
            throw std::invalid_argument( "u values must start at 0 and end at 1" );

        const SweepInput sweep_input = [&sweep]() {
            const auto& m = sweep.mesh;
            std::vector<bool> zero_bcs( m.points.size(), false );
            std::vector<bool> one_bcs( m.points.size(), false );
            for( const VertexId::Type& vid : sweep.source ) zero_bcs.at( vid ) = true;
            for( const VertexId::Type& vid : sweep.target ) one_bcs.at( vid ) = true;
            return SweepInput( m, zero_bcs, one_bcs );
        }();

        using namespace topology;
        using namespace util;

        const size_t degree = 1;
        const basis::KnotVector kv_u( concatenate( {0}, concatenate( u_values, {1} ) ), 1e-9 );
        const basis::KnotVector kv_st = basis::integerKnotsWithNElems( n_elems_st, degree );

        const basis::TPSplineSpace ss = basis::buildBSpline( { kv_st, kv_st, kv_u }, { degree, degree, degree } );
        const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( ss.source() );
        const auto& source_cmap = source_ss.basisComplex().parametricAtlas().cmap();

        const param::ParentDomain pd_3d = param::cubeDomain( 3 );

        const auto to_unit_disk = [&source_cmap]( const topology::Face& f, const Eigen::Vector2d& pt ) -> Eigen::Vector2d {
            const auto [d1, d2, _] = source_cmap.unflatten( f.dart() );
            const Eigen::Vector2d square_coords(
                ( (double)d1.id() + pt( 0 ) ) / (double)cellCount( source_cmap.sourceCMap(), 1 ),
                ( (double)d2.id() + pt( 1 ) ) / (double)cellCount( source_cmap.lineCMap(), 1 ) );

            const double a = 2 * square_coords( 0 ) - 1;
            const double b = 2 * square_coords( 1 ) - 1;
            if( a == 0 and b == 0 ) return Eigen::Vector2d( a, b );
            const double asq = a*a;
            const double bsq = b*b;
            const double scaling = sqrt( asq + bsq - asq * bsq ) / sqrt( asq + bsq );
            return scaling * Eigen::Vector2d( a, b );
        };

        std::map<topology::Vertex, Eigen::Vector2d> tutte_points;
        iterateCellsWhile( source_ss.basisComplex().parametricAtlas().cmap(), 0, [&]( const topology::Vertex& v ) {
            const Eigen::Vector2d pt = source_ss.basisComplex().parametricAtlas().parentPoint( v ).mPoint;
            const Eigen::Vector2d circle_pt = to_unit_disk( topology::Face( v.dart() ), pt );
            tutte_points.emplace( lowestDartId( source_cmap, v ), circle_pt );

            return true;
        } );

        Eigen::MatrixXd fit_cpts;

        reparam::levelSetBasedTracing( sweep_input, u_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
            fit_cpts = reparam::fitLinearMeshToLeaves( ss, leaves, [&]( const topology::Vertex& v ){
                const auto unflattened_cell = unflattenCell( ss.basisComplex().parametricAtlas().cmap(), topology::Vertex( v ) );
                const topology::Vertex v2d(
                    lowestDartId( source_ss.basisComplex().parametricAtlas().cmap(),
                                unflattened_cell.first.value() ) );
                const size_t leaf_ii = unflattened_cell.second.has_value() ? unflattened_cell.second.value().dart().id() : u_values.size() - 1;

                return std::pair<Eigen::Vector2d, size_t>( tutte_points.at( v2d ), leaf_ii );
            }, std::nullopt );
        } );

        if( debug )
        {
            double min_val = std::numeric_limits<double>::infinity();
            double max_val = -std::numeric_limits<double>::infinity();

            eval::SplineSpaceEvaluator evaler( ss, 1 );

            iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
                evaler.localizeElement( vol );
                iterateAdjacentCellsOfRestrictedCell( ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
                    evaler.localizePoint( ss.basisComplex().parametricAtlas().parentPoint( v ) );
                    const double val = evaler.evaluateJacobian( fit_cpts.transpose() ).determinant();
                    min_val = std::min( min_val, val );
                    max_val = std::max( max_val, val );
                    return true;
                } );
                return true;
            } );

            std::cout << "Min jacobian: " << min_val << std::endl;
            std::cout << "Max jacobian: " << max_val << std::endl;

            io::outputBezierMeshToVTK( ss, fit_cpts, "hex_mesh.vtu" );
            io::outputPartialBezierMeshToVTK( ss, fit_cpts, "bad_fit_cells.vtu", [&]( const auto& callback ) {
                iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
                    evaler.localizeElement( vol );
                    double min_jac = std::numeric_limits<double>::infinity();
                    iterateAdjacentCellsOfRestrictedCell( ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
                        evaler.localizePoint( ss.basisComplex().parametricAtlas().parentPoint( v ) );
                        const double val = evaler.evaluateJacobian( fit_cpts.transpose() ).determinant();
                        min_jac = std::min( min_jac, val );
                        return true;
                    } );
                    if( min_jac < 0 )
                    {
                        callback( vol );
                    }
                    return true;
                } );
            } );
        }

        HexMesh out;

        for( Eigen::Index i = 0; i < fit_cpts.rows(); i++ )
        {
            out.points.push_back( fit_cpts.row( i ) );
        }

        iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
            const auto conn = ss.connectivity( vol );
            out.hexes.push_back({});
            std::transform( conn.begin(), conn.end(), out.hexes.back().begin(), []( const auto& c ) {
                return c.id();
            } );
            return true;
        } );

        return out;
    }

    basis::MultiPatchSplineSpace fivePatchSplineSpace( const size_t degree, const basis::KnotVector& kv_u, const size_t n_elems_st )
    {
        using namespace topology;
        const basis::KnotVector kv_st = basis::integerKnotsWithNElems( n_elems_st, degree );

        const std::shared_ptr<const basis::TPSplineSpace> TP_ss = std::make_shared<const basis::TPSplineSpace>(
            basis::buildBSpline( { kv_st, kv_st, kv_u }, { degree, degree, degree } ) );

        const std::array<Dart, 6> connection_darts{ 1,
                                                    19,
                                                    ( 4 * n_elems_st - 3 ) * 6 + 1,
                                                    ( 4 * ( n_elems_st - 1 ) * n_elems_st + 2 ) * 6 + 1,
                                                    ( 4 * ( n_elems_st - 1 ) * n_elems_st + 3 ) * 6 + 1,
                                                    ( n_elems_st * n_elems_st * 4 - 2 ) * 6 + 1 };

        return basis::buildMultiPatchSplineSpace(
            std::vector<std::shared_ptr<const basis::TPSplineSpace>>( 5, TP_ss ),
            {
                { { 0, connection_darts.at( 2 ) }, { 1, connection_darts.at( 1 ) } },
                { { 0, connection_darts.at( 5 ) }, { 2, connection_darts.at( 1 ) } },
                { { 0, connection_darts.at( 4 ) }, { 3, connection_darts.at( 1 ) } },
                { { 0, connection_darts.at( 0 ) }, { 4, connection_darts.at( 1 ) } },
                { { 1, connection_darts.at( 3 ) }, { 2, connection_darts.at( 0 ) } },
                { { 2, connection_darts.at( 3 ) }, { 3, connection_darts.at( 0 ) } },
                { { 3, connection_darts.at( 3 ) }, { 4, connection_darts.at( 0 ) } },
                { { 4, connection_darts.at( 3 ) }, { 1, connection_darts.at( 0 ) } },
            } );
    }

    HexMesh fitFivePatchHexMeshToSweep( const api::Sweep& sweep, const size_t n_elems_st, const std::vector<double>& u_values, const bool debug )
    {
        if( sweep.source.size() == 0 or sweep.target.size() == 0 )
            throw std::invalid_argument( "No source or target surface specified" );
        if( sweep.mesh.points.size() == 0 or sweep.mesh.simplices.size() == 0 )
            throw std::invalid_argument( "No tet mesh supplied" );
        if( u_values.size() < 2 )
            throw std::invalid_argument( "Insufficient u coordinates specified; please include at least two values." );
        if( n_elems_st < 1 )
            throw std::invalid_argument( "Invalid number of elements in the s and t directions specified" );
        if( u_values.front() != 0.0 or u_values.back() != 1.0 )
            throw std::invalid_argument( "u values must start at 0 and end at 1" );

        const SweepInput sweep_input = [&sweep]() {
            const auto& m = sweep.mesh;
            std::vector<bool> zero_bcs( m.points.size(), false );
            std::vector<bool> one_bcs( m.points.size(), false );
            for( const VertexId::Type& vid : sweep.source ) zero_bcs.at( vid ) = true;
            for( const VertexId::Type& vid : sweep.target ) one_bcs.at( vid ) = true;
            return SweepInput( m, zero_bcs, one_bcs );
        }();

        using namespace topology;

        const size_t degree = 1;
        const basis::KnotVector kv_u( util::concatenate( {0}, util::concatenate( u_values, {1} ) ), 1e-9 );

        const basis::MultiPatchSplineSpace mp_ss = fivePatchSplineSpace( degree, kv_u, n_elems_st );

        const std::shared_ptr<const basis::TPSplineSpace> TP_ss = mp_ss.subSpaces().at( 0 );
        const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( TP_ss->source() );

        const param::ParentDomain pd_3d = param::cubeDomain( 3 );

        const auto to_unit_disk =
            [&cmap = source_ss.basisComplex().parametricAtlas().cmap()](
                const size_t patch_ii, const topology::Face& f, const Eigen::Vector2d& pt ) -> Eigen::Vector2d {
            constexpr double a = std::numbers::sqrt2 / 4;
            using std::numbers::pi;
            const auto [d1, d2, _] = cmap.unflatten( f.dart() );
            const Eigen::Vector2d square_coords(
                ( (double)d1.id() + pt( 0 ) ) / (double)cellCount( cmap.sourceCMap(), 1 ),
                ( (double)d2.id() + pt( 1 ) ) / (double)cellCount( cmap.lineCMap(), 1 ) );
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
        };

        std::map<std::pair<size_t, topology::Vertex>, Eigen::Vector2d> tutte_points;
        iterateCellsWhile( source_ss.basisComplex().parametricAtlas().cmap(), 0, [&]( const topology::Vertex& v ) {
            const Eigen::Vector2d pt = source_ss.basisComplex().parametricAtlas().parentPoint( v ).mPoint;

            for( size_t patch_ii = 0; patch_ii < 5; patch_ii++ )
            {
                const Eigen::Vector2d circle_pt = to_unit_disk( patch_ii, topology::Face( v.dart() ), pt );
                tutte_points.emplace( std::pair( patch_ii, lowestDartId( source_ss.basisComplex().parametricAtlas().cmap(), v ) ), circle_pt );
            }

            return true;
        } );

        Eigen::MatrixXd fit_cpts;

        reparam::levelSetBasedTracing( sweep_input, u_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
            fit_cpts = reparam::fitLinearMeshToLeaves( mp_ss, leaves, [&]( const topology::Vertex& v ){
                const auto [patch_id, d_patch] = mp_ss.basisComplex().parametricAtlas().cmap().toLocalDart( v.dart() );
                const auto unflattened_cell = unflattenCell( TP_ss->basisComplex().parametricAtlas().cmap(), topology::Vertex( d_patch ) );
                const topology::Vertex v2d(
                    lowestDartId( source_ss.basisComplex().parametricAtlas().cmap(),
                                unflattened_cell.first.value() ) );
                const size_t leaf_ii = unflattened_cell.second.has_value() ? unflattened_cell.second.value().dart().id() : u_values.size() - 1;

                return std::pair<Eigen::Vector2d, size_t>( tutte_points.at( { patch_id, v2d } ), leaf_ii );
            }, std::nullopt );
        } );

        if( debug )
        {
            double min_val = std::numeric_limits<double>::infinity();
            double max_val = -std::numeric_limits<double>::infinity();

            eval::SplineSpaceEvaluator evaler( mp_ss, 1 );

            iterateCellsWhile( mp_ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
                evaler.localizeElement( vol );
                iterateAdjacentCellsOfRestrictedCell( mp_ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
                    evaler.localizePoint( mp_ss.basisComplex().parametricAtlas().parentPoint( v ) );
                    const double val = evaler.evaluateJacobian( fit_cpts.transpose() ).determinant();
                    min_val = std::min( min_val, val );
                    max_val = std::max( max_val, val );
                    return true;
                } );
                return true;
            } );

            std::cout << "Min jacobian: " << min_val << std::endl;
            std::cout << "Max jacobian: " << max_val << std::endl;

            io::outputBezierMeshToVTK( mp_ss, fit_cpts, "hex_mesh.vtu" );
            io::outputPartialBezierMeshToVTK( mp_ss, fit_cpts, "bad_fit_cells.vtu", [&]( const auto& callback ) {
                iterateCellsWhile( mp_ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
                    evaler.localizeElement( vol );
                    double min_jac = std::numeric_limits<double>::infinity();
                    iterateAdjacentCellsOfRestrictedCell( mp_ss.basisComplex().parametricAtlas().cmap(), vol, 2, 0, [&]( const topology::Vertex& v ) {
                        evaler.localizePoint( mp_ss.basisComplex().parametricAtlas().parentPoint( v ) );
                        const double val = evaler.evaluateJacobian( fit_cpts.transpose() ).determinant();
                        min_jac = std::min( min_jac, val );
                        return true;
                    } );
                    if( min_jac < 0 )
                    {
                        callback( vol );
                    }
                    return true;
                } );
            } );
        }

        HexMesh out;

        for( Eigen::Index i = 0; i < fit_cpts.rows(); i++ )
        {
            out.points.push_back( fit_cpts.row( i ) );
        }

        iterateCellsWhile( mp_ss.basisComplex().parametricAtlas().cmap(), 3, [&]( const topology::Volume& vol ) {
            const auto conn = mp_ss.connectivity( vol );
            out.hexes.push_back({});
            std::transform( conn.begin(), conn.end(), out.hexes.back().begin(), []( const auto& c ) {
                return c.id();
            } );
            return true;
        } );

        return out;
    }
}