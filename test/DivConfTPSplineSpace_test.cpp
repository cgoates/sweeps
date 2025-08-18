#include <catch2/catch_test_macros.hpp>
#include <BSplineSpace1d.hpp>
#include <BasisComplex1d.hpp>
#include <ParametricAtlas1d.hpp>
#include <CombinatorialMap1d.hpp>
#include <KnotVector.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Logging.hpp>
#include <CommonUtils.hpp>
#include <TPSplineSpace.hpp>
#include <VTKOutput.hpp>
#include <IndexOperations.hpp>
#include <VectorConformingBasisComplex.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <VectorConformingHierarchicalTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <VTKOutput.hpp>
#include <iomanip>

using namespace basis;
using namespace param;
using namespace topology;

static constexpr bool OUTPUT_TO_VTK = false;

std::string zeroPadNum( const size_t width, const size_t num )
{
    std::stringstream ss;
    ss << std::setw( width ) << std::setfill( '0' ) << num;
    return ss.str();
}

TEST_CASE( "Simple Div Conf TP Spline space" )
{
    // Input info
    const KnotVector kv1( { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 }, 1e-10 );
    // const KnotVector kv1( {0,0,0,0,1,2,2,2,2}, 1e-10 );
    const KnotVector kv2( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    const TPSplineSpace ss_2d = buildBSpline( {kv1, kv2},  {degree1, degree2} );

    const auto dcbc = std::make_shared<const VectorConformingBasisComplex>( ss_2d.basisComplexPtr() );
    const VectorConformingTPSplineSpace dcss( dcbc, ss_2d );

    const auto& param_2d = dcss.basisComplex().parametricAtlas();
    const auto& cmap_2d = param_2d.cmap();

    Eigen::MatrixX2d geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
    geom.col( 1 ) += sin( 2 * geom.col( 0 ).array() ).matrix();

    if( OUTPUT_TO_VTK ) io::outputBezierMeshToVTK( ss_2d, geom, "bez_test.vtu" );

    eval::SplineSpaceEvaluator primal_evals( ss_2d, 2 );
    eval::SplineSpaceEvaluator vec_evals( dcss, 1 );

    const SmallVector<double, 4> points{ 0.0, 0.3333333, 0.6666666, 1.0 };

    SimplicialComplex output_points;

    std::vector<Eigen::MatrixX3d> vecs( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_2d, 2 ) * points.size() * points.size(), 3 ) );
    std::vector<Eigen::MatrixX3d> vec_dvec( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_2d, 2 ) * points.size() * points.size(), 3 ) );
    size_t i = 0;
    iterateCellsWhile( cmap_2d, 2, [&]( const topology::Face& f ) {
        primal_evals.localizeElement( f );
        vec_evals.localizeElement( f );
        const ParentDomain pd = param_2d.parentDomain( f );
        util::iterateTensorProduct( {points.size(), points.size()}, [&]( const util::IndexVec& indices ){
            const Eigen::Vector2d pt( points.at( indices.at( 0 ) ), points.at( indices.at( 1 ) ) );
            const ParentPoint ppt( pd, pt, {false, false, false, false} );

            primal_evals.localizePoint( ppt );
            vec_evals.localizePoint( ppt );

            const Eigen::VectorXd spatial_point = primal_evals.evaluateManifold( geom.transpose() );
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedHDivBasis( vec_evals, primal_evals, geom.transpose() );//primal_evals.evaluatePiola( geom.transpose() ) * vec_evals.evaluateBasis().transpose();
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedHDivFirstDerivatives( vec_evals, primal_evals, geom.transpose() );

            output_points.simplices.push_back( { output_points.points.size() } );
            output_points.points.push_back( (Eigen::Vector3d() << spatial_point, 0 ).finished() );

            const auto conn = dcss.connectivity( f );

            for( size_t j = 0; j < conn.size(); j++ )
            {
                vecs.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vecs.row( j );
                vec_dvec.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vec_derivs.row( j ).reshaped( 2, 2 ) * param_vecs.row( j ).transpose();
            }
            i++;
        } );
        
        return true;
    } );

    if( OUTPUT_TO_VTK )
    {
        for( size_t func_ii = 0; func_ii < vecs.size(); func_ii++ )
        {
            io::VTKOutputObject out( output_points );
            out.addVertexField( "vecs", vecs.at( func_ii ) );
            out.addVertexField( "derivatives", vec_dvec.at( func_ii ) );
            io::outputSimplicialFieldToVTK( out, "vec_points" + zeroPadNum( 3, func_ii ) + ".vtu" );
        }
    }

    // FIXME: Add asserts!
}

TEST_CASE( "Simple Hierarchical Div Conf TP Spline space" )
{
    // Input info
    // const size_t degree1 = 3;
    // const size_t degree2 = 2;
    // const KnotVector kv1( { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 }, 1e-10 );
    // const KnotVector kv2( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );
    const size_t degree1 = 4;
    const size_t degree2 = 4;
    const KnotVector kv1( { 0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4, 4 }, 1e-10 );
    const KnotVector kv2( { 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2 }, 1e-10 );
    const std::shared_ptr<const TPSplineSpace> level1 = [&](){
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2},  {degree1, degree2} ) );
    }();
    const std::shared_ptr<const TPSplineSpace> level2 = [&](){
        const KnotVector kv1_2 = nAdicRefine( kv1, 3 );
        const KnotVector kv2_2 = nAdicRefine( kv2, 3 );
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1_2, kv2_2},  {degree1, degree2} ) );
    }();
    const std::shared_ptr<const TPSplineSpace> level3 = [&]() {
        const KnotVector kv1_3 = nAdicRefine( kv1, 6 );
        const KnotVector kv2_3 = nAdicRefine( kv2, 6 );
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1_3, kv2_3},  {degree1, degree2} ) );
    }();

    // const HierarchicalTPSplineSpace primal = buildHierarchicalSplineSpace( { level1, level2, level3 }, {
    //     { Face( 0 ), Face( 12 ), Face( 16 ), Face( 20 ), Face( 24 ), Face( 28 ) },
    //     { Face( 8 ), Face( 20 ), Face( 40 ), Face( 44 ), Face( 48 ), Face( 52 ) },
    //     { Face( 24 ), Face( 28 ), Face( 32 ), Face( 36 ), Face( 88 ), Face( 92 ), Face( 96 ), Face( 100 ) }
    // } );

    const HierarchicalTPSplineSpace primal = buildHierarchicalSplineSpace( { level1, level2, level3 }, {
        { Face( 0 ), Face( 12 ), Face( 16 ), Face( 20 ), Face( 24 ), Face( 28 ) },
        { Face( 12 ), Face( 16 ), Face( 28 ), Face( 32 ), Face( 60 ), Face( 64 ), Face( 68 ), Face( 72 ), Face( 76 ), Face( 80 ), Face( 108 ), Face( 112 ), Face( 116 ), Face( 120 ), Face( 124 ), Face( 128 ) },
        { Face( 40 ), Face( 44 ), Face( 48 ), Face( 52 ), Face( 136 ), Face( 140 ), Face( 144 ), Face( 148 ) }
    } );

    const auto dcbc = std::make_shared<const VectorConformingBasisComplex>( primal.basisComplexPtr() );
    const VectorConformingHierarchicalTPSplineSpace dcss( dcbc, primal );

    const auto& param_2d = dcss.basisComplex().parametricAtlas();
    const auto& cmap_2d = param_2d.cmap();

    Eigen::MatrixX2d geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ) } );
    geom.col( 1 ) += sin( 2 * geom.col( 0 ).array() ).matrix();

    // Transfer control points to hierarchical basis
    geom = ( prolongationOperator( primal ) * geom ).eval();

    if( OUTPUT_TO_VTK ) io::outputBezierMeshToVTK( primal, geom, "bez_test.vtu" );

    eval::SplineSpaceEvaluator primal_evals( primal, 2 );
    eval::SplineSpaceEvaluator vec_evals( dcss, 1 );

    const SmallVector<double, 4> points{ 0.0, 0.3333333, 0.6666666, 1.0 };

    SimplicialComplex output_points;

    std::vector<Eigen::MatrixX3d> vecs( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_2d, 2 ) * points.size() * points.size(), 3 ) );
    std::vector<Eigen::MatrixX3d> vec_dvec( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_2d, 2 ) * points.size() * points.size(), 3 ) );
    size_t i = 0;
    iterateCellsWhile( cmap_2d, 2, [&]( const topology::Face& f ) {
    //const topology::Face f( topology::Dart( 0 ) );
        primal_evals.localizeElement( f );
        vec_evals.localizeElement( f );
        const ParentDomain pd = param_2d.parentDomain( f );
        util::iterateTensorProduct( {points.size(), points.size()}, [&]( const util::IndexVec& indices ){
            const Eigen::Vector2d pt( points.at( indices.at( 0 ) ), points.at( indices.at( 1 ) ) );
            const ParentPoint ppt( pd, pt, {false, false, false, false} );

            primal_evals.localizePoint( ppt );
            vec_evals.localizePoint( ppt );

            const Eigen::VectorXd spatial_point = primal_evals.evaluateManifold( geom.transpose() );
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedHDivBasis( vec_evals, primal_evals, geom.transpose() );//primal_evals.evaluatePiola( geom.transpose() ) * vec_evals.evaluateBasis().transpose();
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedHDivFirstDerivatives( vec_evals, primal_evals, geom.transpose() );

            output_points.simplices.push_back( { output_points.points.size() } );
            output_points.points.push_back( (Eigen::Vector3d() << spatial_point, 0 ).finished() );

            const auto conn = dcss.connectivity( f );

            for( size_t j = 0; j < conn.size(); j++ )
            {
                vecs.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vecs.row( j );
                vec_dvec.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vec_derivs.row( j ).reshaped( 2, 2 ) * param_vecs.row( j ).normalized().transpose();
            }
            i++;
        } );
        
        return true;
    } );

    if( OUTPUT_TO_VTK )
    {
        for( size_t func_ii = 0; func_ii < vecs.size(); func_ii++ )
        {
            io::VTKOutputObject out( output_points );
            out.addVertexField( "vecs", vecs.at( func_ii ) );
            out.addVertexField( "derivatives", vec_dvec.at( func_ii ) );
            io::outputSimplicialFieldToVTK( out, "vec_points" + zeroPadNum( 3, func_ii ) + ".vtu" );
        }
    }

    // FIXME: Add asserts!
}

TEST_CASE( "Simple 3d Div Conf TP Spline space" )
{
    // Input info
    const KnotVector kv1( { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 }, 1e-10 );
    // const KnotVector kv1( {0,0,0,0,1,2,2,2,2}, 1e-10 );
    const KnotVector kv2( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    const TPSplineSpace ss_3d = buildBSpline( {kv1, kv2, kv2},  {degree1, degree2, degree2} );

    const auto dcbc = std::make_shared<const VectorConformingBasisComplex>( ss_3d.basisComplexPtr() );
    const VectorConformingTPSplineSpace dcss( dcbc, ss_3d );

    const auto& param_3d = dcss.basisComplex().parametricAtlas();
    const auto& cmap_3d = param_3d.cmap();

    Eigen::MatrixX3d geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv2, degree2 ) } );
    geom.col( 1 ) += sin( 2 * geom.col( 0 ).array() ).matrix();

    if( OUTPUT_TO_VTK ) io::outputBezierMeshToVTK( ss_3d, geom, "bez_test.vtu" );

    eval::SplineSpaceEvaluator primal_evals( ss_3d, 2 );
    eval::SplineSpaceEvaluator vec_evals( dcss, 1 );

    const SmallVector<double, 4> points{ 0.0, 0.3333333, 0.6666666, 1.0 };

    SimplicialComplex output_points;

    std::vector<Eigen::MatrixX3d> vecs( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_3d, 2 ) * points.size() * points.size(), 3 ) );
    std::vector<Eigen::MatrixX3d> vec_dvec( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_3d, 2 ) * points.size() * points.size(), 3 ) );
    size_t i = 0;
    iterateCellsWhile( cmap_3d, 3, [&]( const topology::Volume& f ) {
        primal_evals.localizeElement( f );
        vec_evals.localizeElement( f );
        const ParentDomain pd = param_3d.parentDomain( f );
        util::iterateTensorProduct( {points.size(), points.size(), points.size()}, [&]( const util::IndexVec& indices ){
            const Eigen::Vector3d pt( points.at( indices.at( 0 ) ), points.at( indices.at( 1 ) ), points.at( indices.at( 2 ) ) );
            const ParentPoint ppt( pd, pt, {false, false, false, false, false, false} );

            primal_evals.localizePoint( ppt );
            vec_evals.localizePoint( ppt );

            const Eigen::VectorXd spatial_point = primal_evals.evaluateManifold( geom.transpose() );
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedHDivBasis( vec_evals, primal_evals, geom.transpose() );//primal_evals.evaluatePiola( geom.transpose() ) * vec_evals.evaluateBasis().transpose();
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedHDivFirstDerivatives( vec_evals, primal_evals, geom.transpose() );

            output_points.simplices.push_back( { output_points.points.size() } );
            output_points.points.push_back( spatial_point );

            const auto conn = dcss.connectivity( f );

            for( size_t j = 0; j < conn.size(); j++ )
            {
                vecs.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vecs.row( j );
                vec_dvec.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vec_derivs.row( j ).reshaped( 3, 3 ) * param_vecs.row( j ).transpose();
            }
            i++;
        } );

        return true;
    } );

    if( OUTPUT_TO_VTK )
    {
        for( size_t func_ii = 0; func_ii < vecs.size(); func_ii++ )
        {
            io::VTKOutputObject out( output_points );
            out.addVertexField( "vecs", vecs.at( func_ii ) );
            out.addVertexField( "derivatives", vec_dvec.at( func_ii ) );
            io::outputSimplicialFieldToVTK( out, "vec_points" + zeroPadNum( 3, func_ii ) + ".vtu" );
        }
    }

    // FIXME: Add asserts!
}

TEST_CASE( "Simple 3d Hierarchical Div Conf TP Spline space" )
{
    const size_t degree1 = 4;
    const size_t degree2 = 4;
    const KnotVector kv1( { 0, 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4, 4 }, 1e-10 );
    const KnotVector kv2( { 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2 }, 1e-10 );
    const std::shared_ptr<const TPSplineSpace> level1 = [&](){
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1, kv2, kv2},  {degree1, degree2, degree2} ) );
    }();
    const std::shared_ptr<const TPSplineSpace> level2 = [&](){
        const KnotVector kv1_2 = nAdicRefine( kv1, 3 );
        const KnotVector kv2_2 = nAdicRefine( kv2, 3 );
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1_2, kv2_2, kv2_2},  {degree1, degree2, degree2} ) );
    }();
    const std::shared_ptr<const TPSplineSpace> level3 = [&]() {
        const KnotVector kv1_3 = nAdicRefine( kv1, 6 );
        const KnotVector kv2_3 = nAdicRefine( kv2, 6 );
        return std::make_shared<const TPSplineSpace>( buildBSpline( {kv1_3, kv2_3, kv2_3},  {degree1, degree2, degree2} ) );
    }();

    const std::vector<std::vector<topology::Cell>> leaf_elements = {
        { Volume( 0 ), Volume( 72 ), Volume( 96 ), Volume( 168 ), Volume( 192 ), Volume( 216 ), Volume( 240 ), Volume( 264 ), Volume( 288 ), Volume( 312 ), Volume( 336 ), Volume( 360 ) },
        { Volume(72), Volume(96), Volume(120), Volume(144), Volume(168), Volume(192), Volume(360), Volume(384), Volume(408), Volume(432), Volume(456), Volume(480), Volume(648), Volume(672), Volume(744), Volume(768), Volume(936), Volume(960), Volume(1032), Volume(1056), Volume(1224), Volume(1248), Volume(1272), Volume(1296), Volume(1320), Volume(1344), Volume(1512), Volume(1536), Volume(1560), Volume(1584), Volume(1608), Volume(1632),
          Volume(1800), Volume(1824), Volume(1848), Volume(1872), Volume(1896), Volume(1920), Volume(2088), Volume(2112), Volume(2136), Volume(2160), Volume(2184), Volume(2208), Volume(2376), Volume(2400), Volume(2424), Volume(2448), Volume(2472), Volume(2496), Volume(2664), Volume(2688), Volume(2712), Volume(2736), Volume(2760), Volume(2784), Volume(2952), Volume(2976), Volume(3000), Volume(3024), Volume(3048), Volume(3072), Volume(3240), Volume(3264), Volume(3288), Volume(3312), Volume(3336), Volume(3360),
          Volume(3528), Volume(3552), Volume(3576), Volume(3600), Volume(3624), Volume(3648), Volume(3816), Volume(3840), Volume(3864), Volume(3888), Volume(3912), Volume(3936), Volume(4104), Volume(4128), Volume(4152), Volume(4176), Volume(4200), Volume(4224), Volume(4392), Volume(4416), Volume(4440), Volume(4464), Volume(4488), Volume(4512), Volume(4680), Volume(4704), Volume(4728), Volume(4752), Volume(4776), Volume(4800), Volume(4968), Volume(4992), Volume(5016), Volume(5040), Volume(5064), Volume(5088)
        },
        { Volume(2544), Volume(2568), Volume(2592), Volume(2616), Volume(3120), Volume(3144), Volume(3168), Volume(3192), Volume(3696), Volume(3720), Volume(3744), Volume(3768), Volume(4272), Volume(4296), Volume(4320), Volume(4344), Volume(9456), Volume(9480), Volume(9504), Volume(9528), Volume(10032), Volume(10056), Volume(10080), Volume(10104), Volume(10608), Volume(10632), Volume(10656), Volume(10680), Volume(11184), Volume(11208), Volume(11232), Volume(11256) }
    };
    if( OUTPUT_TO_VTK and false ) // debugging leaf elements
    {
        const Eigen::MatrixXd geom1 = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv2, degree2 ) } );
        const Eigen::MatrixXd geom2 = util::tensorProduct( { grevillePoints( nAdicRefine( kv1, 3 ), degree1 ), grevillePoints( nAdicRefine( kv2, 3 ), degree2 ), grevillePoints( nAdicRefine( kv2, 3 ), degree2 ) } );
        const Eigen::MatrixXd geom3 = util::tensorProduct( { grevillePoints( nAdicRefine( kv1, 6 ), degree1 ), grevillePoints( nAdicRefine( kv2, 6 ), degree2 ), grevillePoints( nAdicRefine( kv2, 6 ), degree2 ) } );
        io::outputPartialBezierMeshToVTK( *level1, geom1, "bez_test_level1.vtu", [&]( std::function<void( const topology::Cell& )> output_cell ) {
            for( const auto& cell : leaf_elements.at( 0 ) )
            {
                output_cell( cell );
            }
        } );
        io::outputPartialBezierMeshToVTK( *level2, geom2, "bez_test_level2.vtu", [&]( std::function<void( const topology::Cell& )> output_cell ) {
            for( const auto& cell : leaf_elements.at( 1 ) )
            {
                output_cell( cell );
            }
        } );
        io::outputPartialBezierMeshToVTK( *level3, geom3, "bez_test_level3.vtu", [&]( std::function<void( const topology::Cell& )> output_cell ) {
            for( const auto& cell : leaf_elements.at( 2 ) )
            {
                output_cell( cell );
            }
        } );
    }
    const HierarchicalTPSplineSpace primal = buildHierarchicalSplineSpace( { level1, level2, level3 }, leaf_elements );

    const auto dcbc = std::make_shared<const VectorConformingBasisComplex>( primal.basisComplexPtr() );
    const VectorConformingHierarchicalTPSplineSpace dcss( dcbc, primal );

    const auto& param_3d = dcss.basisComplex().parametricAtlas();
    const auto& cmap_3d = param_3d.cmap();

    Eigen::MatrixXd geom = util::tensorProduct( { grevillePoints( kv1, degree1 ), grevillePoints( kv2, degree2 ), grevillePoints( kv2, degree2 ) } );
    geom.col( 1 ) += sin( 2 * geom.col( 0 ).array() ).matrix();

    // Transfer control points to hierarchical basis
    geom = ( prolongationOperator( primal ) * geom ).eval();

    if( OUTPUT_TO_VTK ) io::outputBezierMeshToVTK( primal, geom, "bez_test.vtu" );

    geom = geom.transpose().eval();

    eval::SplineSpaceEvaluator primal_evals( primal, 2 );
    eval::SplineSpaceEvaluator vec_evals( dcss, 1 );

    const SmallVector<double, 4> points{ 0.0, 0.3333333, 0.6666666, 1.0 };

    SimplicialComplex output_points;

    std::vector<Eigen::MatrixX3d> vecs( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_3d, 3 ) * points.size() * points.size() * points.size(), 3 ) );
    std::vector<Eigen::MatrixX3d> vec_dvec( dcss.numFunctions(), Eigen::MatrixX3d::Zero( cellCount( cmap_3d, 3 ) * points.size() * points.size() * points.size(), 3 ) );
    size_t i = 0;
    iterateCellsWhile( cmap_3d, 3, [&]( const topology::Volume& f ) {
        primal_evals.localizeElement( f );
        vec_evals.localizeElement( f );
        const ParentDomain pd = param_3d.parentDomain( f );
        util::iterateTensorProduct( {points.size(), points.size(), points.size()}, [&]( const util::IndexVec& indices ){
            const Eigen::Vector3d pt( points.at( indices.at( 0 ) ), points.at( indices.at( 1 ) ), points.at( indices.at( 2 ) ) );
            const ParentPoint ppt( pd, pt, {false, false, false, false, false, false} );

            primal_evals.localizePoint( ppt );
            vec_evals.localizePoint( ppt );

            const Eigen::VectorXd spatial_point = primal_evals.evaluateManifold( geom );
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedHDivBasis( vec_evals, primal_evals, geom );
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedHDivFirstDerivatives( vec_evals, primal_evals, geom );

            output_points.simplices.push_back( { output_points.points.size() } );
            output_points.points.push_back( spatial_point );

            const auto conn = dcss.connectivity( f );

            for( size_t j = 0; j < conn.size(); j++ )
            {
                vecs.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vecs.row( j );
                vec_dvec.at( conn.at( j ) ).row( i ).head( spatial_vecs.cols() ) = spatial_vec_derivs.row( j ).reshaped( 3, 3 ) * param_vecs.row( j ).transpose();
            }
            i++;
        } );

        return true;
    } );

    if( OUTPUT_TO_VTK )
    {
        for( size_t func_ii = 0; func_ii < vecs.size(); func_ii++ )
        {
            io::VTKOutputObject out( output_points );
            out.addVertexField( "vecs", vecs.at( func_ii ) );
            out.addVertexField( "derivatives", vec_dvec.at( func_ii ) );
            io::outputSimplicialFieldToVTK( out, "vec_points" + zeroPadNum( 3, func_ii ) + ".vtu" );
        }
    }

    // FIXME: Add asserts!
}