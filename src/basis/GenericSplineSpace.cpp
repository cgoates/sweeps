#include <GenericSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParametricAtlas.hpp>

namespace basis
{
    GenericSplineSpace::GenericSplineSpace( const BasisComplex& bc,
                                            const std::map<size_t, Eigen::MatrixXd>& ex_ops,
                                            const std::map<size_t, std::vector<FunctionId>>& connectivity )
        : mBasisComplex( bc ), mExtractionOps( ex_ops ), mConnectivity( connectivity )
    {}

    const BasisComplex& GenericSplineSpace::basisComplex() const
    {
        return mBasisComplex;
    }

    Eigen::MatrixXd GenericSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Extraction operators only supported for elements" );
        const auto indexing = indexingOrError( mBasisComplex.parametricAtlas().cmap(), c.dim() );
        return mExtractionOps.at( indexing( c ) );
    }

    std::vector<FunctionId> GenericSplineSpace::connectivity( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Extraction operators only supported for elements" );
        const auto indexing = indexingOrError( mBasisComplex.parametricAtlas().cmap(), c.dim() );
        return mConnectivity.at( indexing( c ) );
    }
}