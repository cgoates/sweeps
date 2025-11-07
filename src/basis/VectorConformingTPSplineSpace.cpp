#include <VectorConformingTPSplineSpace.hpp>
#include <KnotVector.hpp>

namespace basis
{
    SmallVector<std::shared_ptr<const BSplineSpace1d>, 3>
       getReducedDegree1dBases( const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3>& primal_1d_bases )
    {
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> reduced_1d_bases;

        for( const std::shared_ptr<const BSplineSpace1d>& ss_1d : primal_1d_bases )
        {
            const auto reduced_bc = std::make_shared<const BasisComplex1d>( reduceDegree( ss_1d->basisComplex() ) );
            reduced_1d_bases.emplace_back(
                std::make_shared<const BSplineSpace1d>( reduced_bc, reducedOrder( ss_1d->knotVector() ) ) );
        }
        return reduced_1d_bases;
    }

    SmallVector<SmallVector<size_t, 3>, 3> getIndexPatterns( const ConformingType conforming_type, const size_t dim )
    {
        SmallVector<SmallVector<size_t, 3>, 3> index_patterns;

        constexpr size_t PRIMAL = 0;;
        constexpr size_t REDUCED = 1;

        if( dim == 2 )
        {
            index_patterns = ( conforming_type == ConformingType::Divergence )
                ? SmallVector<SmallVector<size_t, 3>, 3>{ { { PRIMAL, REDUCED }, { REDUCED, PRIMAL } } }
                : SmallVector<SmallVector<size_t, 3>, 3>{ { { REDUCED, PRIMAL }, { PRIMAL, REDUCED } } };
        }
        else if( dim == 3 )
        {
            index_patterns =
                ( conforming_type == ConformingType::Divergence )
                    ? SmallVector<SmallVector<size_t, 3>, 3>{
                          { { PRIMAL, REDUCED, REDUCED }, { REDUCED, PRIMAL, REDUCED }, { REDUCED, REDUCED, PRIMAL } } }
                    : SmallVector<SmallVector<size_t, 3>, 3>{
                          { { REDUCED, PRIMAL, PRIMAL }, { PRIMAL, REDUCED, PRIMAL }, { PRIMAL, PRIMAL, REDUCED } } };
        }
        return index_patterns;
    }

    VectorConformingTPSplineSpace::VectorConformingTPSplineSpace( const std::shared_ptr<const VectorConformingBasisComplex>& bc,
                                                                  const TPSplineSpace& primal_basis )
        : mBasisComplex( bc ), mConformingType( bc->conformingType() )
    {
        const size_t dim = primal_basis.basisComplex().parametricAtlas().cmap().dim();

        const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> primal_1d_bases =
            tensorProductComponentSplines( primal_basis );
        if( primal_1d_bases.size() == 0 )
            throw std::runtime_error( "Cannot build VectorConformingTPSplineSpace except over B-spline patch" );

        mReducedDegree1dBases = getReducedDegree1dBases( primal_1d_bases );

        const auto basis_1d_ptr = [&](){
            SmallVector<std::array<std::shared_ptr<const BSplineSpace1d>, 2>, 3> component_1d_bases;
            for( size_t i = 0; i < dim; i++ )
            {
                component_1d_bases.push_back( { primal_1d_bases.at( i ), mReducedDegree1dBases.at( i ) } );
            }
            return [component_1d_bases, index_pattern = getIndexPatterns( mConformingType, dim )]( const size_t vector_component, const size_t tp_index ) {
                return component_1d_bases.at( tp_index ).at( index_pattern.at( vector_component ).at( tp_index ) );
            };
        }();

        if( dim == 2 )
        {
            for( size_t vector_component = 0; vector_component < dim; vector_component++ )
            {
                const auto scalar_tp_bc =
                    std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                      basis_1d_ptr( vector_component, 0 )->basisComplexPtr(),
                                                      basis_1d_ptr( vector_component, 1 )->basisComplexPtr() );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    scalar_tp_bc,
                    basis_1d_ptr( vector_component, 0 ),
                    basis_1d_ptr( vector_component, 1 ) ) );
            }
        }
        else if( dim == 3 )
        {
            const std::shared_ptr<const param::TPParametricAtlas>& source_param =
                dynamic_cast<const TPBasisComplex&>( primal_basis.source().basisComplex() ).parametricAtlasPtr();

            for( size_t vector_component = 0; vector_component < dim; vector_component++ )
            {
                // 2d source
                const auto source_bc =
                    std::make_shared<TPBasisComplex>( source_param,
                                                      basis_1d_ptr( vector_component, 0 )->basisComplexPtr(),
                                                      basis_1d_ptr( vector_component, 1 )->basisComplexPtr() );
                const std::shared_ptr<const TPSplineSpace> source_basis = std::make_shared<TPSplineSpace>(
                    source_bc, basis_1d_ptr( vector_component, 0 ), basis_1d_ptr( vector_component, 1 ) );

                // 3d vector component
                const auto scalar_tp_bc =
                    std::make_shared<TPBasisComplex>( primal_basis.basisComplex().parametricAtlasPtr(),
                                                      source_bc,
                                                      basis_1d_ptr( vector_component, 2 )->basisComplexPtr() );
                mScalarTPBases.push_back( std::make_shared<const TPSplineSpace>(
                    scalar_tp_bc, source_basis, basis_1d_ptr( vector_component, 2 ) ) );
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

        Eigen::MatrixXd result = Eigen::MatrixXd::Zero( rows, cols );
        for( size_t i = 0, row_offset = 0, col_offset = 0; i < scalar_ops.size(); i++ )
        {
            const Eigen::MatrixXd& op = scalar_ops.at( i );
            result.block( row_offset, col_offset, op.rows(), op.cols() ) = op;
            row_offset += op.rows();
            col_offset += op.cols();
        }

        return result;
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