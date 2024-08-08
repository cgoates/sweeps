#include <NavierStokesDiscretization.hpp>
#include <KnotVector.hpp>

namespace api
{
    basis::TPSplineSpace buildH1( const basis::KnotVector& kv_s,
                                  const basis::KnotVector& kv_t,
                                  const size_t degree_s,
                                  const size_t degree_t )
    {
        const auto cmap_1d_s = std::make_shared<const topology::CombinatorialMap1d>( numElements( kv_s ) );
        const auto cmap_1d_t = std::make_shared<const topology::CombinatorialMap1d>( numElements( kv_t ) );
        const auto param_1d_s =
            std::make_shared<const param::ParametricAtlas1d>( cmap_1d_s, parametricLengths( kv_s ) );
        const auto param_1d_t =
            std::make_shared<const param::ParametricAtlas1d>( cmap_1d_t, parametricLengths( kv_t ) );
        const auto bc_1d_s = std::make_shared<const basis::BasisComplex1d>( param_1d_s, degree_s );
        const auto bc_1d_t = std::make_shared<const basis::BasisComplex1d>( param_1d_t, degree_t );
        const auto ss_1d_s = std::make_shared<const basis::BSplineSpace1d>( bc_1d_s, kv_s );
        const auto ss_1d_t = std::make_shared<const basis::BSplineSpace1d>( bc_1d_t, kv_t );
        const auto cmap_2d = std::make_shared<const topology::TPCombinatorialMap>( cmap_1d_s, cmap_1d_t );
        const auto param_2d = std::make_shared<const param::TPParametricAtlas>( cmap_2d, param_1d_s, param_1d_t );
        const auto H1_bc = std::make_shared<const basis::TPBasisComplex>( param_2d, bc_1d_s, bc_1d_t );
        return basis::TPSplineSpace( H1_bc, ss_1d_s, ss_1d_t );
    }

    basis::DivConfTPSplineSpace buildHDIV( const basis::TPSplineSpace& H1 )
    {
        const auto HDIV_bc = std::make_shared<basis::DivConfBasisComplex>( H1.basisComplexPtr() );
        return basis::DivConfTPSplineSpace( HDIV_bc, H1 );
    }

    basis::TPSplineSpace buildL2( const std::shared_ptr<const param::TPParametricAtlas>& param,
                                  const basis::DivConfTPSplineSpace& HDIV )
    {
        const auto L2_bc =
            std::make_shared<basis::TPBasisComplex>( param,
                                                     HDIV.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
                                                     HDIV.reducedDegree1dBases().at( 1 )->basisComplexPtr() );
        return basis::TPSplineSpace( L2_bc, HDIV.reducedDegree1dBases().at( 0 ), HDIV.reducedDegree1dBases().at( 1 ) );
    }

    NavierStokesDiscretization::NavierStokesDiscretization( const basis::KnotVector& kv_s,
                                                            const basis::KnotVector& kv_t,
                                                            const size_t degree_s,
                                                            const size_t degree_t,
                                                            const Eigen::MatrixX2d& cpts )
        : H1_ss( buildH1( kv_s, kv_t, degree_s, degree_t ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap(), { topology::Dart( 0 ) } ),
          cpts( cpts ),
          H1( H1_ss, 2 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}
}