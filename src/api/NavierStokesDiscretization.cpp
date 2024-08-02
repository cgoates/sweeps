#include <NavierStokesDiscretization.hpp>
#include <KnotVector.hpp>

namespace api
{
    NavierStokesDiscretization::NavierStokesDiscretization( const basis::KnotVector& kv_s,
                                                            const basis::KnotVector& kv_t,
                                                            const size_t degree_s,
                                                            const size_t degree_t,
                                                            const Eigen::MatrixX2d& cpts )
        : cmap_1d_s( numElements( kv_s ) ),
          cmap_1d_t( numElements( kv_t ) ),
          param_1d_s( cmap_1d_s, parametricLengths( kv_s ) ),
          param_1d_t( cmap_1d_t, parametricLengths( kv_t ) ),
          bc_1d_s( param_1d_s, degree_s ),
          bc_1d_t( param_1d_t, degree_t ),
          ss_1d_s( bc_1d_s, kv_s ),
          ss_1d_t( bc_1d_t, kv_t ),
          cmap_2d( cmap_1d_s, cmap_1d_t ),
          cmap_bdry( cmap_2d, { topology::Dart( 0 ) } ),
          param_2d( cmap_2d, param_1d_s, param_1d_t ),
          H1_bc( param_2d, bc_1d_s, bc_1d_t ),
          H1_ss( H1_bc, ss_1d_s, ss_1d_t ),
          HDIV_bc( H1_bc ),
          HDIV_ss( HDIV_bc, H1_ss ),
          L2_bc( param_2d,
                 HDIV_ss.reducedDegree1dBases().at( 0 ).basisComplex(),
                 HDIV_ss.reducedDegree1dBases().at( 1 ).basisComplex() ),
          L2_ss( L2_bc, HDIV_ss.reducedDegree1dBases().at( 0 ), HDIV_ss.reducedDegree1dBases().at( 1 ) ),
          cpts( cpts ),
          H1( H1_ss, 2 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}
}