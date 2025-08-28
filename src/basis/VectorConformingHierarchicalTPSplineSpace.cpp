#include <VectorConformingHierarchicalTPSplineSpace.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <KnotVector.hpp>

namespace basis
{
    VectorConformingHierarchicalTPSplineSpace::VectorConformingHierarchicalTPSplineSpace(
        const std::shared_ptr<const VectorConformingBasisComplex>& bc, const HierarchicalTPSplineSpace& primal_basis )
        : mBasisComplex( bc )
    {
        // Construct d HierarchicalTPSplineSpaces for the d vector components by constructing a series of
        // VectorConformingTPSplineSpaces, and then pulling out the underlying scalar bases.

        if( &bc->parametricAtlas() != &primal_basis.basisComplex().parametricAtlas() )
            throw std::invalid_argument(
                "The basis complex and primal basis given to VectorConformingHierarchicalTPSplineSpace must have the same parametric atlas" );

        const size_t num_levels = primal_basis.basisComplex().parametricAtlas().cmap().numLevels();
        const size_t dim = primal_basis.basisComplex().parametricAtlas().cmap().dim();

        SmallVector<std::vector<std::shared_ptr<const TPBasisComplex>>, 3> scalar_level_bcs(
            dim, std::vector<std::shared_ptr<const TPBasisComplex>>() );
        SmallVector<std::vector<std::shared_ptr<const TPSplineSpace>>, 3> scalar_level_bases(
            dim, std::vector<std::shared_ptr<const TPSplineSpace>>() );
        for( auto& v : scalar_level_bases ) v.reserve( num_levels );
        for( auto& v : scalar_level_bcs ) v.reserve( num_levels );

        const auto& primal_refinement_levels = primal_basis.refinementLevels();

        for( size_t level = 0; level < num_levels; level++ )
        {
            const auto level_bc = std::make_shared<const VectorConformingBasisComplex>( primal_refinement_levels.at( level )->basisComplexPtr(), bc->conformingType() );
            const VectorConformingTPSplineSpace level_ss( level_bc, *primal_refinement_levels.at( level ), bc->conformingType() );
            const SmallVector<std::shared_ptr<const TPSplineSpace>,3>& scalar = level_ss.scalarTPBases();
            if( scalar.size() != scalar_level_bases.size() ) throw std::runtime_error( "Scalar size is not the same" );

            for( size_t i = 0; i < scalar.size(); i++ )
            {
                scalar_level_bcs.at( i ).push_back( scalar.at( i )->basisComplexPtr() );
                scalar_level_bases.at( i ).push_back( scalar.at( i ) );
            }
        }

        const auto& param = primal_basis.basisComplex().parametricAtlasPtr();
        for( size_t i = 0; i < scalar_level_bases.size(); i++ )
        {
            const auto scalar_bc = std::make_shared<const HierarchicalTPBasisComplex>( param, scalar_level_bcs.at( i ) );
            mScalarTPBases.push_back( std::make_shared<const HierarchicalTPSplineSpace>( scalar_bc, scalar_level_bases.at( i ) ) );
        }
    }

    const VectorConformingBasisComplex& VectorConformingHierarchicalTPSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd VectorConformingHierarchicalTPSplineSpace::extractionOperator( const topology::Cell& c ) const
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

    std::vector<FunctionId> VectorConformingHierarchicalTPSplineSpace::connectivity( const topology::Cell& c ) const
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

    size_t VectorConformingHierarchicalTPSplineSpace::numFunctions() const
    {
        return std::accumulate( mScalarTPBases.begin(), mScalarTPBases.end(), 0, [&]( const size_t accum, const auto& ss ) {
            return accum + ss->numFunctions();
        } );
    }
}