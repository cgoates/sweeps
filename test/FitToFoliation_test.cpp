#include <catch2/catch_test_macros.hpp>
#include <SimplicialComplexTestCases.hpp>
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
#include <Foliation.hpp>

using util::linspace;
using util::concatenate;
using namespace reparam;

Eigen::Vector2d toUnitSquare( const topology::TPCombinatorialMap& cmap, const topology::Face& f, const Eigen::Vector2d& pt )
{
    // FIXME: Assumes unit parametric lengths
    const auto [d1, d2, _] = cmap.unflatten( f.dart() );
    return Eigen::Vector2d( ( (double)d1.id() + pt( 0 ) ) / (double)cellCount( cmap.sourceCMap(), 1 ), ( (double)d2.id() + pt( 1 ) ) / (double)cellCount( cmap.lineCMap(), 1 ) );
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

Eigen::MatrixXd fitToLeaves(
    const SweepInput& sweep_input,
    eval::SplineSpaceEvaluator& evaler,
    const std::vector<double>& level_set_values,
    const size_t n_points,
    const std::function<void(
        const size_t,
        const std::function<Eigen::Vector3d( const topology::Cell&, const param::ParentPoint&, const Eigen::Vector2d& )>& )>&
        leaf_point_iterator )
{
    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    // If not covering the whole domain, add source and target surface
    std::vector<double> level_set_values2;
    const bool add_front = not util::equals( level_set_values.front(), 0.0, 1e-9 );
    const bool add_back = not util::equals( level_set_values.back(), 1.0, 1e-9 );
    if( add_front ) level_set_values2.push_back( 0 );
    level_set_values2.insert( level_set_values2.end(), level_set_values.begin(), level_set_values.end() );
    if( add_back ) level_set_values2.push_back( 1.0 );

    Eigen::MatrixXd fit_cpts;

    reparam::levelSetBasedTracing( sweep_input, level_set_values2, [&]( const std::vector<FoliationLeaf>& leaves2 ) {
        const std::vector<FoliationLeaf> leaves( add_front ? std::next( leaves2.begin() ) : leaves2.begin(), add_back ? std::prev( leaves2.end() ) : leaves2.end() );
        std::cout << "About to fit\n";

        fit_cpts = fitting::leastSquaresFitting(
            evaler,
            n_points,
            3,
            [&]( const std::function<void( const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd& )>&
                    add_least_squares_point ) {
                for( size_t leaf_ii = 0; leaf_ii < leaves.size(); leaf_ii++ )
                {
                    leaf_point_iterator( leaf_ii, [&]( const topology::Cell& cell, const param::ParentPoint& vol_ppt, const Eigen::Vector2d& circle_pt ) {
                        const auto param_pt = leaves.at( leaf_ii ).tutte_mapping->maybeInverse( circle_pt );
                        CHECK( param_pt.has_value() );
                        if( not param_pt.has_value() )
                        {
                            std::cerr << "NO VALUE" << std::endl;
                            pauseDebugger();
                        }
                        const Eigen::Vector3d field_pt = leaves.at( leaf_ii ).space_mapping->evaluate( param_pt.value().first, param_pt.value().second );

                        add_least_squares_point( cell, vol_ppt, field_pt );
                        return field_pt;
                    } );
                }
            } );
    } );

    return fit_cpts;
}

void fitToPringlesSinglePatch( const SweepInput& sweep_input,
                               const std::vector<double>& level_set_values,
                               const std::string& output_prefix,
                               const size_t degree,
                               const size_t n_elems_u,
                               const size_t n_elems_st )
{
    const basis::KnotVector kv1 = basis::integerKnotsWithNElems( n_elems_st, degree );
    const basis::KnotVector kv2 = basis::integerKnotsWithNElems( n_elems_u, degree );

    const basis::TPSplineSpace vol_ss = basis::buildBSpline( {kv1, kv1, kv2}, {degree, degree, degree} );
    const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( vol_ss.source() );
    const topology::TPCombinatorialMap& vol_cmap = vol_ss.basisComplex().parametricAtlas().cmap();

    const std::vector<std::pair<topology::Cell, param::ParentPoint>> ppt_u =
        parentPointsOfParamPoints( level_set_values, vol_ss.line().basisComplex().parametricAtlas(), 1e-9 );

    const size_t n_points = cellCount( source_ss.basisComplex().parametricAtlas().cmap(), 2 ) * 4 * level_set_values.size();

    const SmallVector<double, 2> source_points{ 0.1, 0.9 };

    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    std::cout << "About to fit\n";

    SimplicialComplex fitting_points;

    eval::SplineSpaceEvaluator evaler( vol_ss, 0 );

    const Eigen::MatrixXd fit_cpts = fitToLeaves( sweep_input, evaler, level_set_values, n_points, [&]( const size_t leaf_ii, const auto& point_callback ) {
        iterateCellsWhile(
            source_ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const topology::Cell& f ) {
                util::iterateTensorProduct(
                    { source_points.size(), source_points.size() }, [&]( const util::IndexVec& indices ) {
                        const Eigen::Vector3d pt( source_points.at( indices.at( 0 ) ),
                                                source_points.at( indices.at( 1 ) ),
                                                ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
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

                        const Eigen::Vector3d field_pt = point_callback( cell, vol_ppt, circle_pt );

                        fitting_points.simplices.emplace_back( fitting_points.points.size() );
                        fitting_points.points.push_back( field_pt );
                    } );
                return true;
            } );
    } );

    io::VTKOutputObject fitting_points_output( fitting_points );
    io::outputSimplicialFieldToVTK( fitting_points_output, output_prefix + "fitting_points.vtu" );

    std::cout << "Finished fit\n";

    std::cout << "Vol cells: " << cellCount( vol_cmap, 3 ) << std::endl;

    io::outputBezierMeshToVTK( vol_ss,
                                fit_cpts,
                                "fit_to_" + output_prefix + ".vtu" );
}

void fitToPringles5Patch( const SweepInput& sweep_input,
                          const std::vector<double>& level_set_values,
                          const std::string& output_prefix,
                          const size_t degree,
                          const basis::KnotVector& kv_u,
                          const size_t n_elems_st )
{
    using namespace topology;

    const basis::KnotVector kv_st = basis::integerKnotsWithNElems( n_elems_st, degree );

    const std::shared_ptr<const basis::TPSplineSpace> TP_ss = std::make_shared<const basis::TPSplineSpace>(
        basis::buildBSpline( { kv_st, kv_st, kv_u }, { degree, degree, degree } ) );
    const auto& vol_cmap = TP_ss->basisComplex().parametricAtlas().cmap();
    const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( TP_ss->source() );
    const std::vector<std::pair<topology::Cell, param::ParentPoint>> ppt_u =
        parentPointsOfParamPoints( level_set_values, TP_ss->line().basisComplex().parametricAtlas(), 1e-9 );

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

    
    const SmallVector<double, 2> source_points{ 0.1, 0.9 };
    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    SimplicialComplex fitting_points;
    eval::SplineSpaceEvaluator evaler( mp_ss, 0 );

    const size_t n_bdry_pts = 8 * level_set_values.size();

    const size_t n_points = cellCount( source_ss.basisComplex().parametricAtlas().cmap(), 2 ) * 4 * 5 * level_set_values.size() + n_bdry_pts;

    Eigen::MatrixXd fit_cpts = fitToLeaves( sweep_input, evaler, level_set_values, n_points, [&]( const size_t leaf_ii, const auto& point_callback ) {
        for( size_t patch_ii = 0; patch_ii < 5; patch_ii++ )
        {
            iterateCellsWhile(
                source_ss.basisComplex().parametricAtlas().cmap(), 2, [&]( const topology::Cell& f ) {
                    util::iterateTensorProduct(
                        { source_points.size(), source_points.size() }, [&]( const util::IndexVec& indices ) {
                            const Eigen::Vector3d pt( source_points.at( indices.at( 0 ) ),
                                                    source_points.at( indices.at( 1 ) ),
                                                    ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
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
                            
                            const Eigen::Vector3d field_pt = point_callback( cell, vol_ppt, circle_pt );
                            fitting_points.simplices.emplace_back( fitting_points.points.size() );
                            fitting_points.points.push_back( field_pt );
                        } );
                    return true;
                } );


            if( patch_ii > 1 )
            {
                for( const topology::Edge& e : { topology::Edge( 5 ), topology::Edge( 13 ) } )
                {
                    const param::ParentPoint surf_ppt = param::pointOnBoundary( param::cubeDomain( 2 ), param::parentDomainBoundary( source_ss.basisComplex().parametricAtlas(), e ) );
                    const Eigen::Vector3d pt( surf_ppt.mPoint( 0 ), surf_ppt.mPoint( 1 ), ppt_u.at( leaf_ii ).second.mPoint( 0 ) );
                    const param::ParentPoint vol_ppt(
                        pd_3d,
                        pt,
                        { surf_ppt.mBaryCoordIsZero.at( 0 ),
                        surf_ppt.mBaryCoordIsZero.at( 1 ),
                        surf_ppt.mBaryCoordIsZero.at( 2 ),
                        surf_ppt.mBaryCoordIsZero.at( 3 ),
                        ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 0 ),
                        ppt_u.at( leaf_ii ).second.mBaryCoordIsZero.at( 1 ) } );

                    const topology::Volume patch_cell(
                        vol_cmap.flatten( e.dart(),
                                        ppt_u.at( leaf_ii ).first.dart(),
                                        topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );

                    const topology::Volume cell( mp_ss.basisComplex().parametricAtlas().cmap().toGlobalDart( patch_ii, patch_cell.dart() ) );

                    const Eigen::Vector2d circle_pt = multiPatchToUnitDisk(
                        patch_ii,
                        toUnitSquare(
                            source_ss.basisComplex().parametricAtlas().cmap(), topology::Face( e.dart() ), pt.head( 2 ) ) );

                    const Eigen::Vector3d field_pt = point_callback( cell, vol_ppt, circle_pt );
                    fitting_points.simplices.emplace_back( fitting_points.points.size() );
                    fitting_points.points.push_back( field_pt );
                }
            }
        }
    } );

    io::VTKOutputObject fitting_points_output( fitting_points );
    io::outputSimplicialFieldToVTK( fitting_points_output, output_prefix + "fitting_points.vtu" );

    io::outputBezierMeshToVTK( mp_ss, fit_cpts, "fit_to_" + output_prefix + "_multi_patch.vtu" );
}

void fitToPringles5Patch( const SweepInput& sweep_input,
                          const std::vector<double>& level_set_values,
                          const std::string& output_prefix,
                          const size_t degree,
                          const size_t& n_elems_u,
                          const size_t n_elems_st )
{
    const basis::KnotVector kv_u = basis::integerKnotsWithNElems( n_elems_u, degree );
    return fitToPringles5Patch( sweep_input, level_set_values, output_prefix, degree, kv_u, n_elems_st );
}

void linearMeshPringles5Patch( const SweepInput& sweep_input,
                               const std::string& output_prefix,
                               const std::vector<double>& level_set_values,
                               const size_t n_elems_st )
{
    using namespace topology;

    const size_t degree = 1;
    const basis::KnotVector kv_u( concatenate( {0}, concatenate( level_set_values, {1} ) ), 1e-9 );
    const basis::KnotVector kv_st = basis::integerKnotsWithNElems( n_elems_st, degree );

    const std::shared_ptr<const basis::TPSplineSpace> TP_ss = std::make_shared<const basis::TPSplineSpace>(
        basis::buildBSpline( { kv_st, kv_st, kv_u }, { degree, degree, degree } ) );
    const basis::TPSplineSpace& source_ss = static_cast<const basis::TPSplineSpace&>( TP_ss->source() );

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

    const param::ParentDomain pd_3d = param::cubeDomain( 3 );

    std::map<std::pair<size_t, topology::Vertex>, Eigen::Vector2d> tutte_points;
    iterateCellsWhile( source_ss.basisComplex().parametricAtlas().cmap(), 0, [&]( const topology::Vertex& v ) {
        const Eigen::Vector2d pt = source_ss.basisComplex().parametricAtlas().parentPoint( v ).mPoint;

        const Eigen::Vector2d sq_pt = toUnitSquare( source_ss.basisComplex().parametricAtlas().cmap(), topology::Face( v.dart() ), pt );

        for( size_t patch_ii = 0; patch_ii < 5; patch_ii++ )
        {
            const Eigen::Vector2d circle_pt = multiPatchToUnitDisk( patch_ii, sq_pt );
            tutte_points.emplace( std::pair( patch_ii, lowestDartId( source_ss.basisComplex().parametricAtlas().cmap(), v ) ), circle_pt );
        }

        return true;
    } );

    Eigen::MatrixXd fit_cpts;

    reparam::levelSetBasedTracing( sweep_input, level_set_values, [&]( const std::vector<FoliationLeaf>& leaves ) {
        fit_cpts = reparam::fitLinearMeshToLeaves( mp_ss, leaves, [&]( const topology::Vertex& v ){
            const auto [patch_id, d_patch] = mp_ss.basisComplex().parametricAtlas().cmap().toLocalDart( v.dart() );
            const auto unflattened_cell = unflattenCell( TP_ss->basisComplex().parametricAtlas().cmap(), topology::Vertex( d_patch ) );
            const topology::Vertex v2d(
                lowestDartId( source_ss.basisComplex().parametricAtlas().cmap(),
                              unflattened_cell.first.value() ) );
            const size_t leaf_ii = unflattened_cell.second.has_value() ? unflattened_cell.second.value().dart().id() : level_set_values.size() - 1;

            return std::pair<Eigen::Vector2d, size_t>( tutte_points.at( { patch_id, v2d } ), leaf_ii );
        }, std::nullopt );
    } );

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

    io::outputBezierMeshToVTK( mp_ss, fit_cpts, "fit_to_" + output_prefix + "_multi_patch.vtu" );
    io::outputPartialBezierMeshToVTK( mp_ss, fit_cpts, "bad_fit_cells_" + output_prefix + ".vtu", [&]( const auto& callback ) {
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
    const std::string output_prefix = "ventricle";

    // fitToPringlesSinglePatch( sweep_input, level_set_values, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, output_prefix, 2, 20, 2 );
}

TEST_CASE( "Linear mesh on left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );

    // level set every 0.02 until 0.75, then ~30 levels between 0.75 and 0.85, then every 0.2 until 1.0.
    // const std::vector<double> level_set_values =
    //     concatenate(
    //         concatenate( linspace( 0, 0.78, 40 ), linspace( 0.8, 0.819, 20 ) ),
    //         concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 31 ) ), linspace( 0.86, 1.0, 13 ) ) );
    const std::vector<double> level_set_values =
        concatenate(
            concatenate( linspace( 0, 0.78, 20 ), linspace( 0.8, 0.819, 10 ) ),
            concatenate( concatenate( linspace( 0.81925, 0.81975, 3 ), linspace( 0.82, 0.85, 15 ) ), linspace( 0.86, 1.0, 7 ) ) );
    // const std::vector<double> level_set_values = linspace( 0, 1.0, 80 );
    const std::string output_prefix = "ventricle";

    // fitToPringlesSinglePatch( sweep_input, level_set_values, output_prefix, 2, 20, 4 );
    linearMeshPringles5Patch( sweep_input, output_prefix, level_set_values, 2 );
}

TEST_CASE( "Level set parameterization of part of left ventricle" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/left_ventricle.inp", "Surface3", "Surface2" );
    // const std::vector<double> level_set_values = linspace( 0.75, 0.85, 80 );
    const std::vector<double> level_set_values =
        concatenate( concatenate( linspace( 0.75, 0.8, 4 ), linspace( 0.8, 0.818, 14 ) ),
                     concatenate( linspace( 0.81825, 0.82275, 19 ), linspace( 0.823, 0.85, 31 ) ) );
    const std::string output_prefix = "part_ventricle";

    const basis::KnotVector kv_u( {
        0, 0, 0, 6, 12, 12.7, 13, 13.25, 13.5, 13.75, 14, 14.25, 14.5, 14.75, 15, 15.4, 16, 17, 18, 19, 20, 20, 20 
    }, 1e-9 );
    fitToPringles5Patch( sweep_input, level_set_values, output_prefix, 2, kv_u, 2 );
}

TEST_CASE( "Level set parameterization of hook" )
{
    const SweepInput sweep_input = io::loadINPFile( SRC_HOME "/test/data/hook.inp", "Surface12", "Surface10" );

    const std::string output_prefix = "hook";
    const std::vector<double> level_set_values = linspace( 0, 1.0, 100 );

    const basis::KnotVector kv_u( concatenate( { 0, 0, 0, 0.015, 0.015, 0.16, 0.16515 },
                                               concatenate( linspace( 0.173, 0.73, 15 ), concatenate( linspace( 0.75, 1.0, 5 ), { 1.0, 1.0 } ) ) ),
                                  1e-9 );

    // fitToPringlesSinglePatch( sweep_input, level_set_values, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, output_prefix, 2, kv_u, 2 );
}

TEST_CASE( "Level set parameterization of femur" )
{
    const SweepInput sweep_input = [&](){
        SweepInput sweep_input =
            io::loadINPFile( SRC_HOME "/test/data/femur.inp", "Surface5", "Surface2" );
        for( size_t i = 0; i < sweep_input.mesh.points.size(); i++ )
            sweep_input.one_bcs.at( i ) = sweep_input.mesh.points.at( i )( 1 ) > -30;
        return sweep_input;
    }();

    const std::string output_prefix = "femur";
    const std::vector<double> level_set_values = linspace( 0, 1.0, 100 );

    // const basis::KnotVector kv_u = basis::integerKnotsWithNElems( 20, 2 );
    const basis::KnotVector kv_u( concatenate( {0, 0}, concatenate( linspace( 0, 16, 9 ), concatenate( linspace( 16.5, 20, 15 ), {20, 20} ) ) ), 1e-9 );

    // fitToPringlesSinglePatch( sweep_input, level_set_values, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, output_prefix, 2, kv_u, 2 );
}

TEST_CASE( "Level set parameterization of spring" )
{
    const SweepInput sweep_input =
            io::loadINPFile( SRC_HOME "/test/data/spring.inp", "Surface2", "Surface3" );

    const std::string output_prefix = "spring";
    const std::vector<double> level_set_values = linspace( 0, 1.0, 200 );

    const basis::KnotVector kv_u = basis::integerKnotsWithNElems( 100, 2 );
    // const basis::KnotVector kv_u( concatenate( {0}, concatenate( linspace( 0, 16, 9 ), concatenate( linspace( 16.5, 20, 15 ), {20} ) ) ), 1e-9 );

    // fitToPringlesSinglePatch( sweep_input, level_set_values, output_prefix, 2, 20, 4 );
    fitToPringles5Patch( sweep_input, level_set_values, output_prefix, 2, kv_u, 2 );
}