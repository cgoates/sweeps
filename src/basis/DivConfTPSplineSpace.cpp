#include <DivConfTPSplineSpace.hpp>
#include <KnotVector.hpp>

namespace basis
{
    DivConfTPSplineSpace::DivConfTPSplineSpace( const DivConfBasisComplex& bc, const TPSplineSpace& primal_basis ) :
        mBasisComplex( bc )
    {
        const size_t dim = primal_basis.basisComplex().parametricAtlas().cmap().dim();
        mReducedDegree1dBasisComplex.reserve( dim );
        mReducedDegree1dBases.reserve( dim );
        mScalarTPBases.reserve( dim );

        std::vector<std::reference_wrapper<const BSplineSpace1d>> primal_1d_bases;
        primal_1d_bases.reserve( dim );

        // Pull out the 1d basis complexes and degree reduce them
        primal_1d_bases.push_back( std::cref( primal_basis.line() ) );
        if( dim == 3 )
        {
            const TPSplineSpace& source_primal = static_cast<const TPSplineSpace&>( primal_basis.source() );
            primal_1d_bases.push_back( source_primal.line() );
            primal_1d_bases.push_back( static_cast<const BSplineSpace1d&>( source_primal.source() ) );
        }
        else
            primal_1d_bases.push_back( static_cast<const BSplineSpace1d&>( primal_basis.source() ) );
        
        for( const BSplineSpace1d& ss_1d : primal_1d_bases )
        {
            mReducedDegree1dBasisComplex.push_back( reduceDegree( ss_1d.basisComplex() ) );
            mReducedDegree1dBases.emplace_back( mReducedDegree1dBasisComplex.back(), reducedOrder( ss_1d.knotVector() ) );
        }

        if( dim == 2 )
        {
            mScalarTPBasisComplexes.reserve( dim );
            mScalarTPBasisComplexes.emplace_back( primal_basis.basisComplex().parametricAtlas(),
                                                  primal_1d_bases.at( 0 ).get().basisComplex(),
                                                  mReducedDegree1dBasisComplex.at( 1 ) );
            mScalarTPBases.emplace_back( mScalarTPBasisComplexes.back(), primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
            mScalarTPBasisComplexes.emplace_back( primal_basis.basisComplex().parametricAtlas(),
                                                  mReducedDegree1dBasisComplex.at( 0 ),
                                                  primal_1d_bases.at( 1 ).get().basisComplex() );
            mScalarTPBases.emplace_back( mScalarTPBasisComplexes.back(), mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) );
        }
        else if( dim == 3 )
        {
            mScalarTPBasisComplexes.reserve( 2 * dim );
            m2dSourceTPBases.reserve( dim );

            const param::TPParametricAtlas& source_param = static_cast<const param::TPParametricAtlas&>( primal_basis.source().basisComplex().parametricAtlas() );

            // 2d source of component 0
            mScalarTPBasisComplexes.emplace_back( source_param,
                                                  primal_1d_bases.at( 0 ).get().basisComplex(),
                                                  mReducedDegree1dBasisComplex.at( 1 ) );
            m2dSourceTPBases.emplace_back( mScalarTPBasisComplexes.back(), primal_1d_bases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
            // component 0
            mScalarTPBasisComplexes.emplace_back( primal_basis.basisComplex().parametricAtlas(), mScalarTPBasisComplexes.back(), mReducedDegree1dBasisComplex.at( 2 ) );
            mScalarTPBases.emplace_back( mScalarTPBasisComplexes.back(), m2dSourceTPBases.back(), mReducedDegree1dBases.at( 2 ) );

            // 2d source of component 1
            mScalarTPBasisComplexes.emplace_back( source_param,
                                                  mReducedDegree1dBasisComplex.at( 0 ),
                                                  primal_1d_bases.at( 1 ).get().basisComplex() );
            m2dSourceTPBases.emplace_back( mScalarTPBasisComplexes.back(), mReducedDegree1dBases.at( 0 ), primal_1d_bases.at( 1 ) );
            // component 1
            mScalarTPBasisComplexes.emplace_back( primal_basis.basisComplex().parametricAtlas(), mScalarTPBasisComplexes.back(), mReducedDegree1dBasisComplex.at( 2 ) );
            mScalarTPBases.emplace_back( mScalarTPBasisComplexes.back(), m2dSourceTPBases.back(), mReducedDegree1dBases.at( 2 ) );

            // 2d source of component 1
            mScalarTPBasisComplexes.emplace_back( source_param,
                                                  mReducedDegree1dBasisComplex.at( 0 ),
                                                  mReducedDegree1dBasisComplex.at( 1 ) );
            m2dSourceTPBases.emplace_back( mScalarTPBasisComplexes.back(), mReducedDegree1dBases.at( 0 ), mReducedDegree1dBases.at( 1 ) );
            // component 1
            mScalarTPBasisComplexes.emplace_back( primal_basis.basisComplex().parametricAtlas(), mScalarTPBasisComplexes.back(), primal_1d_bases.at( 2 ).get().basisComplex() );
            mScalarTPBases.emplace_back( mScalarTPBasisComplexes.back(), m2dSourceTPBases.back(), primal_1d_bases.at( 2 ) );
        }
    }

    const DivConfBasisComplex& DivConfTPSplineSpace::basisComplex() const
    {
        return mBasisComplex;
    }

    Eigen::MatrixXd DivConfTPSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        std::vector<Eigen::MatrixXd> scalar_ops;
        scalar_ops.reserve( mScalarTPBases.size() );
        std::transform( mScalarTPBases.begin(), mScalarTPBases.end(), std::back_inserter( scalar_ops ), [&]( const TPSplineSpace& scalar_basis ) {
            return scalar_basis.extractionOperator( c );
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
                     Eigen::MatrixXd::Zero( scalar_ops.at( 1 ).rows(), scalar_ops.at( 1 ).cols() ),
                     Eigen::MatrixXd::Zero( scalar_ops.at( 0 ).rows(), scalar_ops.at( 0 ).cols() ),
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
        for( const TPSplineSpace& scalar_basis : mScalarTPBases )
        {
            const std::vector<FunctionId> scalar_connectivity = scalar_basis.connectivity( c );
            connectivity.reserve( connectivity.size() + scalar_connectivity.size() );

            std::transform( scalar_connectivity.begin(),
                            scalar_connectivity.end(),
                            std::back_inserter( connectivity ),
                            [&]( const FunctionId& fid ) { return FunctionId( fid + offset ); } );

            offset += scalar_connectivity.size();
        }

        return connectivity;
    }

    size_t DivConfTPSplineSpace::numFunctions() const
    {
        return std::accumulate( mScalarTPBases.begin(), mScalarTPBases.end(), 0, [&]( const size_t accum, const TPSplineSpace& ss ) {
            return accum + ss.numFunctions();
        } );
    }
}