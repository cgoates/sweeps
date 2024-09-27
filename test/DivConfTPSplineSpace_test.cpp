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
#include <DivConfBasisComplex.hpp>
#include <DivConfTPSplineSpace.hpp>
#include <DivConfHierarchicalTPSplineSpace.hpp>
#include <SplineSpaceEvaluator.hpp>
#include <VTKOutput.hpp>
#include <format>

using namespace basis;
using namespace param;
using namespace topology;

static constexpr bool OUTPUT_TO_VTK = false;

TEST_CASE( "Simple Div Conf TP Spline space" )
{
    // Input info
    const KnotVector kv1( { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 }, 1e-10 );
    // const KnotVector kv1( {0,0,0,0,1,2,2,2,2}, 1e-10 );
    const KnotVector kv2( { 0, 0, 0, 1, 2, 2, 2 }, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    const TPSplineSpace ss_2d = buildBSpline( {kv1, kv2},  {degree1, degree2} );

    const auto dcbc = std::make_shared<const DivConfBasisComplex>( ss_2d.basisComplexPtr() );
    const DivConfTPSplineSpace dcss( dcbc, ss_2d );

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
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedVectorBasis( vec_evals, primal_evals, geom.transpose() );//primal_evals.evaluatePiola( geom.transpose() ) * vec_evals.evaluateBasis().transpose();
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedVectorFirstDerivatives( vec_evals, primal_evals, geom.transpose() );

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
            io::outputSimplicialFieldToVTK( out, "vec_points" + std::format( "{:0>3}", func_ii ) + ".vtu" );
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

    const auto dcbc = std::make_shared<const DivConfBasisComplex>( primal.basisComplexPtr() );
    const DivConfHierarchicalTPSplineSpace dcss( dcbc, primal );

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
            const Eigen::MatrixXd spatial_vecs = eval::piolaTransformedVectorBasis( vec_evals, primal_evals, geom.transpose() );//primal_evals.evaluatePiola( geom.transpose() ) * vec_evals.evaluateBasis().transpose();
            const Eigen::MatrixXd param_vecs = vec_evals.evaluateBasis();
            const Eigen::MatrixXd spatial_vec_derivs = eval::piolaTransformedVectorFirstDerivatives( vec_evals, primal_evals, geom.transpose() );

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
            io::outputSimplicialFieldToVTK( out, "vec_points" + std::format( "{:0>3}", func_ii ) + ".vtu" );
        }
    }

    // FIXME: Add asserts!
}