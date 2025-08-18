#include <VectorConformingTPSplineSpace.hpp>
#include <KnotVector.hpp>

namespace basis
{
    VectorConformingTPSplineSpace::VectorConformingTPSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc,
                                                const TPSplineSpace& primal_basis,
                                                const ConformingType conforming_type )
        : mBasisComplex( bc ), mConformingType( conforming_type )
    {
        const size_t dim = primal_basis.basisComplex().parametricAtlas().cmap().dim();

        const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> primal_1d_bases =
            tensorProductComponentSplines( primal_basis );
        if( primal_1d_bases.size() == 0 )
            throw std::runtime_error( "Cannot build VectorConformingTPSplineSpace except over B-spline patch" );

        SmallVector<std::shared_ptr<const BasisComplex1d>, 3> reduced_1d_bcs;
        std::shared_ptr<const TPBasisComplex> temp_tp_bc;

        for( const std::shared_ptr<const BSplineSpace1d>& ss_1d : primal_1d_bases )
        {
            reduced_1d_bcs.push_back( std::make_shared<const BasisComplex1d>( reduceDegree( ss_1d->basisComplex() ) ) );
            mReducedDegree1dBases.emplace_back(
                std::make_shared<const BSplineSpace1d>( reduced_1d_bcs.back(), reducedOrder( ss_1d->knotVector() ) ) );
        }

        if( dim == 2 )
        {
            if( conforming_type == ConformingType::Divergence )
            {
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               primal_1d_bases.at( 0 )->basisComplexPtr(),
                                                               reduced_1d_bcs.at( 1 ) );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) ) );
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               reduced_1d_bcs.at( 0 ),
                                                               primal_1d_bases.at( 1 )->basisComplexPtr() );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) ) );
            }
            else if( conforming_type == ConformingType::Curl )
            {
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               reduced_1d_bcs.at( 0 ),
                                                               primal_1d_bases.at( 1 )->basisComplexPtr() );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) ) );
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               primal_1d_bases.at( 0 )->basisComplexPtr(),
                                                               reduced_1d_bcs.at( 1 ) );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) ) );
            }
        }
        else if( dim == 3 )
        {
            std::shared_ptr<const TPSplineSpace> temp_tp_2d_basis;

            const std::shared_ptr<const param::TPParametricAtlas>& source_param =
                dynamic_cast<const TPBasisComplex&>( primal_basis.source().basisComplex() ).parametricAtlasPtr();

            if( conforming_type == ConformingType::Divergence )
            {
                // 2d source of component 0
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    source_param, primal_1d_bases.at( 0 )->basisComplexPtr(), reduced_1d_bcs.at( 1 ) );
                temp_tp_2d_basis = std::make_shared<TPSplineSpace>(
                    temp_tp_bc, primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
                // component 0
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    primal_basis.basisComplex().parametricAtlasPtr(), temp_tp_bc, reduced_1d_bcs.at( 2 ) );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, temp_tp_2d_basis, mReducedDegree1dBases.at( 2 ) ) );

                // 2d source of component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    source_param, reduced_1d_bcs.at( 0 ), primal_1d_bases.at( 1 )->basisComplexPtr() );
                temp_tp_2d_basis = std::make_shared<TPSplineSpace>(
                    temp_tp_bc, mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) );
                // component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    primal_basis.basisComplex().parametricAtlasPtr(), temp_tp_bc, reduced_1d_bcs.at( 2 ) );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, temp_tp_2d_basis, mReducedDegree1dBases.at( 2 ) ) );

                // 2d source of component 1
                temp_tp_bc =
                    std::make_shared<TPBasisComplex>( source_param, reduced_1d_bcs.at( 0 ), reduced_1d_bcs.at( 1 ) );
                temp_tp_2d_basis = std::make_shared<TPSplineSpace>(
                    temp_tp_bc, mReducedDegree1dBases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
                // component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               temp_tp_bc,
                                                               primal_1d_bases.at( 2 )->basisComplexPtr() );
                mScalarTPBases.push_back(
                    std::make_shared<const TPSplineSpace>( temp_tp_bc, temp_tp_2d_basis, primal_1d_bases.at( 2 ) ) );
            }
            else if( conforming_type == ConformingType::Curl )
            {
                // 2d source of component 0
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    source_param, reduced_1d_bcs.at( 0 ), primal_1d_bases.at( 1 )->basisComplexPtr() );
                temp_tp_2d_basis = std::make_shared<TPSplineSpace>(
                    temp_tp_bc, mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) );
                // component 0
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               temp_tp_bc,
                                                               primal_1d_bases.at( 2 )->basisComplexPtr() );
                mScalarTPBases.push_back(
                    std::make_shared<const TPSplineSpace>( temp_tp_bc, temp_tp_2d_basis, primal_1d_bases.at( 2 ) ) );

                // 2d source of component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    source_param, primal_1d_bases.at( 0 )->basisComplexPtr(), reduced_1d_bcs.at( 1 ) );
                temp_tp_2d_basis = std::make_shared<TPSplineSpace>(
                    temp_tp_bc, primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
                // component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                               temp_tp_bc,
                                                               primal_1d_bases.at( 2 )->basisComplexPtr() );
                mScalarTPBases.push_back(
                    std::make_shared<const TPSplineSpace>( temp_tp_bc, temp_tp_2d_basis, primal_1d_bases.at( 2 ) ) );

                // 2d source of component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>( source_param,
                                                               primal_1d_bases.at( 0 )->basisComplexPtr(),
                                                               primal_1d_bases.at( 1 )->basisComplexPtr() );
                temp_tp_2d_basis =
                    std::make_shared<TPSplineSpace>( temp_tp_bc, primal_1d_bases.at( 0 ), primal_1d_bases.at( 1 ) );
                // component 1
                temp_tp_bc = std::make_shared<TPBasisComplex>(
                    primal_basis.basisComplex().parametricAtlasPtr(), temp_tp_bc, reduced_1d_bcs.at( 2 ) );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    temp_tp_bc, temp_tp_2d_basis, mReducedDegree1dBases.at( 2 ) ) );
            }
        }
    }

    const VectorConformingBasisComplex& VectorConformingTPSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd VectorConformingTPSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        SmallVector<Eigen::MatrixXd, 3> scalar_ops;
        std::transform( mScalarTPBases.begin(), mScalarTPBases.end(), std::back_inserter( scalar_ops ), [&]( const auto& scalar_basis ) {
            return scalar_basis->extractionOperator( c );
        } );
        const auto [rows, cols] = std::accumulate( scalar_ops.begin(),
                                                   scalar_ops.end(),
                                                   std::pair<Eigen::Index, Eigen::Index>( 0, 0 ),
                                                   [&]( const auto& accum, const Eigen::MatrixXd& op ) {
                                                       return std::pair<Eigen::Index, Eigen::Index>(
                                                           accum.first + op.rows(), accum.second + op.cols() );
                                                   } );

        if( scalar_ops.size() == 2 )
        {
            return ( Eigen::MatrixXd( rows, cols ) << scalar_ops.at( 0 ),
                     Eigen::MatrixXd::Zero( scalar_ops.at( 0 ).rows(), scalar_ops.at( 1 ).cols() ),
                     Eigen::MatrixXd::Zero( scalar_ops.at( 1 ).rows(), scalar_ops.at( 0 ).cols() ),
                     scalar_ops.at( 1 ) )
                .finished();
        }
        else if( scalar_ops.size() == 3 )
        {
            const auto zeros = [&]( const size_t r, const size_t c ) { return Eigen::MatrixXd::Zero( scalar_ops.at( r ).rows(), scalar_ops.at( c ).cols() ); };
            return ( Eigen::MatrixXd( rows, cols ) << scalar_ops.at( 0 ), zeros( 0, 1 ), zeros( 0, 2 ),
                                                      zeros( 1, 0 ), scalar_ops.at( 1 ), zeros( 1, 2 ),
                                                      zeros( 2, 0 ), zeros( 2, 1 ), scalar_ops.at( 2 ) ).finished();
        }
        throw std::runtime_error( "Bad dimension for div conf tp spline space" );
    }

    std::vector<FunctionId> VectorConformingTPSplineSpace::connectivity( const topology::Cell& c ) const
    {
        std::vector<FunctionId> connectivity;

        size_t offset = 0;
        for( const auto& scalar_basis : mScalarTPBases )
        {
            const std::vector<FunctionId> scalar_connectivity = scalar_basis->connectivity( c );
            connectivity.reserve( connectivity.size() + scalar_connectivity.size() );

            std::transform( scalar_connectivity.begin(),
                            scalar_connectivity.end(),
                            std::back_inserter( connectivity ),
                            [&]( const FunctionId& fid ) { return FunctionId( fid + offset ); } );

            offset += scalar_basis->numFunctions();
        }

        return connectivity;
    }

    size_t VectorConformingTPSplineSpace::numFunctions() const
    {
        return std::accumulate( mScalarTPBases.begin(), mScalarTPBases.end(), 0, [&]( const size_t accum, const auto& ss ) {
            return accum + ss->numFunctions();
        } );
    }

    size_t VectorConformingTPSplineSpace::numVectorComponents() const
    {
        return mScalarTPBases.size();
    }

}