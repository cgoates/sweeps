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
#include <SplineSpaceEvaluator.hpp>
#include <VTKOutput.hpp>
#include <format>

using namespace basis;
using namespace param;
using namespace topology;

TEST_CASE( "TP Spline space evaluation" )
{
    // Input info
    const KnotVector kv1( {0,0,0,0,1,2,3,4,4,4,4}, 1e-10 );
    const KnotVector kv2( {0,0,0,1,2,2,2}, 1e-10 );
    const size_t degree1 = 3;
    const size_t degree2 = 2;

    // Built from input info
    const CombinatorialMap1d cmap1( numElements( kv1 ) );
    const CombinatorialMap1d cmap2( numElements( kv2 ) );
    const ParametricAtlas1d param1( cmap1, parametricLengths( kv1 ) );
    const ParametricAtlas1d param2( cmap2, parametricLengths( kv2 ) );
    const BasisComplex1d bc1( param1, degree1 );
    const BasisComplex1d bc2( param2, degree2 );
    const BSplineSpace1d ss1( bc1, kv1 );
    const BSplineSpace1d ss2( bc2, kv2 );

    const TPCombinatorialMap cmap_2d( cmap1, cmap2 );
    const TPParametricAtlas param_2d( cmap_2d, param1, param2 );
    const TPBasisComplex bc_2d( param_2d, bc1, bc2 );
    const TPSplineSpace ss_2d( bc_2d, ss1, ss2 );

    eval::SplineSpaceEvaluator primal_evals( ss_2d, 1 );

    iterateCellsWhile( cmap_2d, 1, [&]( const topology::Edge& e ) {
        const auto maybe_opp_e = phi( cmap_2d, 2, e.dart() );
        if( maybe_opp_e.has_value() )
        {
            const topology::Face f( e.dart() );
            const topology::Face f_opp( maybe_opp_e.value() );
            const ParentPoint ppt = pointOnBoundary( param_2d.parentDomain( f ), parentDomainBoundary( param_2d, e ) );
            primal_evals.localizeElement( f );
            primal_evals.localizePoint( ppt );

            const auto conn = ss_2d.connectivity( f );
            const auto basis = primal_evals.evaluateBasis();
            const auto derv = primal_evals.evaluateFirstDerivatives();

            const ParentPoint ppt_opp = pointOnBoundary( param_2d.parentDomain( f_opp ), parentDomainBoundary( param_2d, topology::Edge( maybe_opp_e.value() ) ) );
            primal_evals.localizeElement( f_opp );
            primal_evals.localizePoint( ppt_opp );

            const auto conn_opp = ss_2d.connectivity( f_opp );
            const auto basis_opp = primal_evals.evaluateBasis();
            const auto derv_opp = primal_evals.evaluateFirstDerivatives();

            for( size_t i = 0; i < conn.size(); i++ )
            {
                const auto it = std::find( conn_opp.begin(), conn_opp.end(), conn.at( i ) );
                if( it != conn_opp.end() )
                {
                    const size_t j = std::distance( conn_opp.begin(), it );

                    CHECK( util::equals( basis.row( i ), basis_opp.row( j ), 1e-9 ) );
                    CHECK( util::equals( derv.row( i ), derv_opp.row( j ), 1e-9 ) );

                    if( not util::equals( derv.row( i ), derv_opp.row( j ), 1e-9 ) )
                    {
                        std::cout << e.dart() << ", " << conn.at( i ) << ":\n";
                        std::cout << basis.row( i ) << " vs " << basis_opp.row( j ) << std::endl;
                        std::cout << derv.row( i ) << " vs " << derv_opp.row( j ) << std::endl;
                        std::cout << "----------\n";
                    }
                }
            }
        }
        return true;
    } );
}