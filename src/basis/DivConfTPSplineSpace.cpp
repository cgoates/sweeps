#include <DivConfTPSplineSpace.hpp>
#include <KnotVector.hpp>

namespace basis
{
    DivConfTPSplineSpace::DivConfTPSplineSpace( const std::shared_ptr<const DivConfBasisComplex>& bc,
                                                const TPSplineSpace& primal_basis )
        : mBasisComplex( bc )
    {
        const size_t dim = primal_basis.basisComplex().parametricAtlas().cmap().dim();

        const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> primal_1d_bases =
            tensorProductComponentSplines( primal_basis );
        if( primal_1d_bases.size() == 0 )
            throw std::runtime_error( "Cannot build DivConfTPSplineSpace except over B-spline patch" );

        SmallVector<std::shared_ptr<const BasisComplex1d>, 3> reduced_1d_bcs;
        SmallVector<std::shared_ptr<const TPBasisComplex>, 3> scalar_tp_bcs;

        for( const std::shared_ptr<const BSplineSpace1d>& ss_1d : primal_1d_bases )
        {
            reduced_1d_bcs.push_back( std::make_shared<const BasisComplex1d>( reduceDegree( ss_1d->basisComplex() ) ) );
            mReducedDegree1dBases.emplace_back(
                std::make_shared<const BSplineSpace1d>( reduced_1d_bcs.back(), reducedOrder( ss_1d->knotVector() ) ) );
        }

        if( dim == 2 )
        {
            scalar_tp_bcs.emplace_back(
                std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                  primal_1d_bases.at( 0 )->basisComplexPtr(),
                                                  reduced_1d_bcs.at( 1 ) ) );
            mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>( scalar_tp_bcs.back(), primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) ) );
            scalar_tp_bcs.emplace_back(
                std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                  reduced_1d_bcs.at( 0 ),
                                                  primal_1d_bases.at( 1 )->basisComplexPtr() ) );
            mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>( scalar_tp_bcs.back(), mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) ) );
        }
        else if( dim == 3 )
        {
            SmallVector<std::shared_ptr<const TPSplineSpace>, 3> tp_2d_bases;

            const std::shared_ptr<const param::TPParametricAtlas>& source_param =
                dynamic_cast<const TPBasisComplex&>( primal_basis.source().basisComplex() ).parametricAtlasPtr();

            // 2d source of component 0
            scalar_tp_bcs.emplace_back( std::make_shared<TPBasisComplex>(
                source_param, primal_1d_bases.at( 0 )->basisComplexPtr(), reduced_1d_bcs.at( 1 ) ) );
            tp_2d_bases.emplace_back( std::make_shared<TPSplineSpace>(
                scalar_tp_bcs.back(), primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) ) );
            // component 0
            scalar_tp_bcs.emplace_back( std::make_shared<TPBasisComplex>(
                primal_basis.basisComplex().parametricAtlasPtr(), scalar_tp_bcs.back(), reduced_1d_bcs.at( 2 ) ) );
            mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>( scalar_tp_bcs.back(), tp_2d_bases.back(), mReducedDegree1dBases.at( 2 ) ) );

            // 2d source of component 1
            scalar_tp_bcs.emplace_back( std::make_shared<TPBasisComplex>(
                source_param, reduced_1d_bcs.at( 0 ), primal_1d_bases.at( 1 )->basisComplexPtr() ) );
            tp_2d_bases.emplace_back( std::make_shared<TPSplineSpace>(
                scalar_tp_bcs.back(), mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) ) );
            // component 1
            scalar_tp_bcs.emplace_back( std::make_shared<TPBasisComplex>(
                primal_basis.basisComplex().parametricAtlasPtr(), scalar_tp_bcs.back(), reduced_1d_bcs.at( 2 ) ) );
            mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>( scalar_tp_bcs.back(), tp_2d_bases.back(), mReducedDegree1dBases.at( 2 ) ) );

            // 2d source of component 1
            scalar_tp_bcs.emplace_back(
                std::make_shared<TPBasisComplex>( source_param, reduced_1d_bcs.at( 0 ), reduced_1d_bcs.at( 1 ) ) );
            tp_2d_bases.emplace_back( std::make_shared<TPSplineSpace>(
                scalar_tp_bcs.back(), mReducedDegree1dBases.at( 0 ), mReducedDegree1dBases.at( 1 ) ) );
            // component 1
            scalar_tp_bcs.emplace_back(
                std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                  scalar_tp_bcs.back(),
                                                  primal_1d_bases.at( 2 )->basisComplexPtr() ) );
            mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>( scalar_tp_bcs.back(), tp_2d_bases.back(), primal_1d_bases.at( 2 ) ) );
        }
    }

    const DivConfBasisComplex& DivConfTPSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd DivConfTPSplineSpace::extractionOperator( const topology::Cell& c ) const
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
            const Eigen::MatrixXd zeros_0 = Eigen::MatrixXd::Zero( scalar_ops.at( 0 ).rows(), scalar_ops.at( 0 ).cols() );
            const Eigen::MatrixXd zeros_1 = Eigen::MatrixXd::Zero( scalar_ops.at( 1 ).rows(), scalar_ops.at( 1 ).cols() );
            const Eigen::MatrixXd zeros_2 = Eigen::MatrixXd::Zero( scalar_ops.at( 2 ).rows(), scalar_ops.at( 2 ).cols() );
            return ( Eigen::MatrixXd( rows, cols ) << scalar_ops.at( 0 ), zeros_1, zeros_2,
                                                      zeros_0, scalar_ops.at( 1 ), zeros_2,
                                                      zeros_0, zeros_1, scalar_ops.at( 2 ) ).finished();
        }
        throw std::runtime_error( "Bad dimension for div conf tp spline space" );
    }

    std::vector<FunctionId> DivConfTPSplineSpace::connectivity( const topology::Cell& c ) const
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

    size_t DivConfTPSplineSpace::numFunctions() const
    {
        return std::accumulate( mScalarTPBases.begin(), mScalarTPBases.end(), 0, [&]( const size_t accum, const auto& ss ) {
            return accum + ss->numFunctions();
        } );
    }
}