#include <FittingUtilities.hpp>
#include <SplineSpace.hpp>
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>
#include <CombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>

namespace fitting
{

Eigen::MatrixXd linearControlPointsFromVertexPositions( const basis::SplineSpace& ss, const VertexPositionsFunc& positions )
{
    Eigen::MatrixXd fit_cpts = Eigen::MatrixXd::Zero( cellCount( ss.basisComplex().parametricAtlas().cmap(), 0 ), 3 );

    iterateCellsWhile( ss.basisComplex().parametricAtlas().cmap(), 0, [&]( const topology::Vertex& v ) {
        const auto conn = ss.connectivity( v );
        if( conn.size() != 1 )
        {
            throw std::runtime_error( "Vertex " + std::to_string( v.dart().id() ) + " has " +
                                      std::to_string( conn.size() ) + " functions, expected 1." );
        }

        fit_cpts.row( conn.front() ) = positions( v );
        return true;
    } );

    return fit_cpts;
}

}