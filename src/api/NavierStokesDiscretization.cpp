#include <NavierStokesDiscretization.hpp>
#include <KnotVector.hpp>

namespace api
{
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

    NavierStokesTPDiscretization::NavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                                                const basis::KnotVector& kv_t,
                                                                const size_t degree_s,
                                                                const size_t degree_t,
                                                                const Eigen::Matrix2Xd& cpts )
        : H1_ss( basis::buildBSpline( { kv_s, kv_t }, { degree_s, degree_t } ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap(), { topology::Dart( 0 ) } ),
          cpts( cpts ),
          H1( H1_ss, 2 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}
}