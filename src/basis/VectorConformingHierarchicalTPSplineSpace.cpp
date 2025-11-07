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
            const VectorConformingTPSplineSpace level_ss( level_bc, *primal_refinement_levels.at( level ) );
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

        Eigen::MatrixXd out = Eigen::MatrixXd::Zero( rows, cols );
        for( size_t i = 0, row_offset = 0, col_offset = 0; i < scalar_ops.size(); i++ )
        {
            const Eigen::MatrixXd& op = scalar_ops.at( i );
            out.block( row_offset, col_offset, op.rows(), op.cols() ) = op;
            row_offset += op.rows();
            col_offset += op.cols();
        }

        return out;
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