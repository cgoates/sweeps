#include <VectorConformingHierarchicalTPSplineSpace.hpp>
#include <VectorConformingTPSplineSpace.hpp>
#include <HierarchicalTPParametricAtlas.hpp>
#include <KnotVector.hpp>
#include <algorithm>

namespace basis
{
    SmallVector<std::vector<std::shared_ptr<const TPSplineSpace>>, 3>
        scalarRefinementLevels( const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& refinement_levels )
    {
        if( refinement_levels.empty() )
            throw std::invalid_argument( "VectorConformingHierarchicalTPSplineSpace requires at least one refinement level." );

        const size_t num_components = refinement_levels.front()->numVectorComponents();
        SmallVector<std::vector<std::shared_ptr<const TPSplineSpace>>, 3> scalar_level_bases(
            num_components, std::vector<std::shared_ptr<const TPSplineSpace>>() );
        for( auto& bases : scalar_level_bases ) bases.reserve( refinement_levels.size() );

        for( const auto& level : refinement_levels )
        {
            if( level->numVectorComponents() != num_components )
                throw std::invalid_argument(
                    "All vector conforming refinement levels must have the same number of components." );

            const SmallVector<std::shared_ptr<const TPSplineSpace>, 3>& scalar = level->scalarTPBases();
            if( scalar.size() != num_components )
                throw std::runtime_error( "Vector conforming level has inconsistent scalar component count." );

            for( size_t component = 0; component < num_components; component++ )
            {
                scalar_level_bases.at( component ).push_back( scalar.at( component ) );
            }
        }

        return scalar_level_bases;
    }

    std::vector<size_t> componentOffsets( const VectorConformingTPSplineSpace& level )
    {
        std::vector<size_t> offsets;
        offsets.reserve( level.scalarTPBases().size() + 1 );
        offsets.push_back( 0 );
        for( const auto& scalar_basis : level.scalarTPBases() )
        {
            offsets.push_back( offsets.back() + scalar_basis->numFunctions() );
        }
        if( offsets.back() != level.numFunctions() )
            throw std::runtime_error( "Vector conforming component offsets do not span the whole level." );
        return offsets;
    }

    std::vector<std::vector<std::vector<FunctionId>>>
        splitActiveVectorFunctionsByComponent(
            const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& refinement_levels,
            const std::vector<std::vector<FunctionId>>& active_funcs )
    {
        if( refinement_levels.size() != active_funcs.size() )
            throw std::invalid_argument(
                "Active vector function lists must have one entry per vector refinement level." );

        const size_t num_components = refinement_levels.front()->numVectorComponents();
        std::vector<std::vector<std::vector<FunctionId>>> scalar_active_funcs(
            num_components, std::vector<std::vector<FunctionId>>( refinement_levels.size() ) );

        for( size_t level_ii = 0; level_ii < refinement_levels.size(); level_ii++ )
        {
            const auto& level = *refinement_levels.at( level_ii );
            const std::vector<size_t> offsets = componentOffsets( level );

            for( const FunctionId& fid : active_funcs.at( level_ii ) )
            {
                if( static_cast<size_t>( fid.id() ) >= level.numFunctions() )
                    throw std::invalid_argument( "Active vector function id is outside its refinement level." );

                const auto upper =
                    std::upper_bound( offsets.begin(), offsets.end(), static_cast<size_t>( fid.id() ) );
                if( upper == offsets.begin() )
                    throw std::runtime_error( "Could not locate component for active vector function." );

                const size_t component = static_cast<size_t>( std::distance( offsets.begin(), upper ) - 1 );
                scalar_active_funcs.at( component ).at( level_ii ).push_back(
                    FunctionId( fid.id() - offsets.at( component ) ) );
            }
        }

        return scalar_active_funcs;
    }

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

    VectorConformingHierarchicalTPSplineSpace::VectorConformingHierarchicalTPSplineSpace(
        const std::shared_ptr<const VectorConformingBasisComplex>& bc,
        const std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>>& refinement_levels,
        const std::vector<std::vector<FunctionId>>& active_funcs )
        : mBasisComplex( bc )
    {
        if( refinement_levels.empty() )
            throw std::invalid_argument( "VectorConformingHierarchicalTPSplineSpace requires refinement levels." );

        const auto primal_hier_bc =
            std::dynamic_pointer_cast<const HierarchicalTPBasisComplex>( bc->primalComplexPtr() );
        if( primal_hier_bc == nullptr )
            throw std::invalid_argument(
                "Explicit active VectorConformingHierarchicalTPSplineSpace requires a HierarchicalTPBasisComplex." );

        if( primal_hier_bc->parametricAtlas().cmap().numLevels() != refinement_levels.size() )
            throw std::invalid_argument(
                "VectorConformingHierarchicalTPSplineSpace basis complex and refinement levels disagree on level count." );

        const auto scalar_level_bases = scalarRefinementLevels( refinement_levels );
        const auto scalar_active_funcs = splitActiveVectorFunctionsByComponent( refinement_levels, active_funcs );

        for( size_t component = 0; component < scalar_level_bases.size(); component++ )
        {
            std::vector<std::shared_ptr<const TPBasisComplex>> scalar_level_bcs;
            scalar_level_bcs.reserve( scalar_level_bases.at( component ).size() );
            for( const auto& level_basis : scalar_level_bases.at( component ) )
            {
                scalar_level_bcs.push_back( level_basis->basisComplexPtr() );
            }

            const auto scalar_bc =
                std::make_shared<const HierarchicalTPBasisComplex>(
                    primal_hier_bc->parametricAtlasPtr(), scalar_level_bcs );
            mScalarTPBases.push_back(
                std::make_shared<const HierarchicalTPSplineSpace>(
                    scalar_bc, scalar_level_bases.at( component ), scalar_active_funcs.at( component ) ) );
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
