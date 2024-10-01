#include <TPSplineSpace.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <IndexOperations.hpp>

namespace basis
{
    TPSplineSpace::TPSplineSpace( const std::shared_ptr<const TPBasisComplex>& bc,
                                  const std::shared_ptr<const SplineSpace>& source,
                                  const std::shared_ptr<const BSplineSpace1d>& line )
        : mBasisComplex( bc ), mSource( source ), mLine( line )
    {
        if( &source->basisComplex() != &mBasisComplex->sourceComplex() or
            &line->basisComplex() != &mBasisComplex->lineComplex() )
            throw std::runtime_error( "Invalid inputs to TPSplineSpace::TPSplineSpace()" );
    }

    const TPBasisComplex& TPSplineSpace::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd TPSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex->parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Bad cell dimension for spline space" );

        const auto dart_pr = mBasisComplex->parametricAtlas().cmap().unflatten( c.dart() );
        const topology::Cell source_c( std::get<0>( dart_pr ), mSource->basisComplex().parametricAtlas().cmap().dim() );
        const topology::Cell line_c( std::get<1>( dart_pr ), mLine->basisComplex().parametricAtlas().cmap().dim() );
        return Eigen::kroneckerProduct( mLine->extractionOperator( line_c ), mSource->extractionOperator( source_c ) );
    }

    FunctionId TPSplineSpace::flatten( const FunctionId& source_fid, const FunctionId& line_fid ) const
    {
        return FunctionId( source_fid + mSource->numFunctions() * line_fid );
    }


    std::vector<FunctionId> TPSplineSpace::connectivity( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex->parametricAtlas().cmap().dim() )
            throw std::runtime_error( "Bad cell dimension for spline space" );

        const auto dart_pr = mBasisComplex->parametricAtlas().cmap().unflatten( c.dart() );
        const std::vector<FunctionId> source_conn = mSource->connectivity(
            topology::Cell( std::get<0>( dart_pr ), mSource->basisComplex().parametricAtlas().cmap().dim() ) );
        const std::vector<FunctionId> line_conn = mLine->connectivity(
            topology::Cell( std::get<1>( dart_pr ), mLine->basisComplex().parametricAtlas().cmap().dim() ) );

        std::vector<FunctionId> out;
        out.reserve( line_conn.size() * source_conn.size() );
        for( const FunctionId& line_fid : line_conn )
        {
            for( const FunctionId& source_fid : source_conn )
            {
                out.push_back( flatten( source_fid, line_fid ) );
            }
        }

        return out;
    }

    size_t TPSplineSpace::numFunctions() const
    {
        return mSource->numFunctions() * mLine->numFunctions();
    }

    SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> tensorProductComponentSplines( const TPSplineSpace& ss )
    {
        const size_t dim = ss.basisComplex().parametricAtlas().cmap().dim();
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> out;

        if( dim == 3 )
        {
            const std::shared_ptr<const TPSplineSpace> source_primal =
                std::dynamic_pointer_cast<const TPSplineSpace>( ss.sourcePtr() );
            if( source_primal.get() == nullptr ) return {};
            out.push_back( std::dynamic_pointer_cast<const BSplineSpace1d>( source_primal->sourcePtr() ) );
            if( out.back().get() == nullptr ) return {};
            out.push_back( source_primal->linePtr() );
        }
        else
        {
            out.push_back( std::dynamic_pointer_cast<const BSplineSpace1d>( ss.sourcePtr() ) );
            if( out.back().get() == nullptr ) return {};
        }
        out.push_back( ss.linePtr() );
        return out;
    }

    util::IndexVec getTPLengths( const basis::TPSplineSpace& tpss )
    {
        const auto component_comps = tensorProductComponentSplines( tpss );
        util::IndexVec out;
        for( const auto& comp : component_comps ) out.push_back( comp->numFunctions() );
        return out;
    };

    TPSplineSpace buildBSpline( const SmallVector<KnotVector, 3>& kvs, const SmallVector<size_t, 3>& degrees )
    {
        SmallVector<std::shared_ptr<const topology::CombinatorialMap1d>, 3> cmap_1ds;
        SmallVector<std::shared_ptr<const param::ParametricAtlas1d>, 3> atlas_1ds;
        SmallVector<std::shared_ptr<const BasisComplex1d>, 3> bc_1ds;
        SmallVector<std::shared_ptr<const BSplineSpace1d>, 3> ss_1ds;

        for( size_t i = 0; i < kvs.size(); i++ )
        {
            const auto& kv = kvs.at( i );
            const size_t deg = degrees.at( i );
            cmap_1ds.push_back( std::make_shared<const topology::CombinatorialMap1d>( numElements( kv ) ) );
            atlas_1ds.push_back( std::make_shared<const param::ParametricAtlas1d>( cmap_1ds.back(), parametricLengths( kv ) ) );
            bc_1ds.push_back( std::make_shared<const BasisComplex1d>( atlas_1ds.back(), deg ) );
            ss_1ds.push_back( std::make_shared<const BSplineSpace1d>( bc_1ds.back(), kv ) );
        }

        auto cmap = std::make_shared<const topology::TPCombinatorialMap>( cmap_1ds.at( 0 ), cmap_1ds.at( 1 ) );
        auto atlas = std::make_shared<const param::TPParametricAtlas>( cmap, atlas_1ds.at( 0 ), atlas_1ds.at( 1 ) );
        auto bc = std::make_shared<const TPBasisComplex>( atlas, bc_1ds.at( 0 ), bc_1ds.at( 1 ) );

        if( kvs.size() == 2 )
        {
            return TPSplineSpace( bc, ss_1ds.at( 0 ), ss_1ds.at( 1 ) );
        }
        else if( kvs.size() == 3 )
        {
            const auto ss = std::make_shared<const TPSplineSpace>( bc, ss_1ds.at( 0 ), ss_1ds.at( 1 ) );

            cmap = std::make_shared<const topology::TPCombinatorialMap>( cmap, cmap_1ds.at( 2 ) );
            atlas = std::make_shared<const param::TPParametricAtlas>( cmap, atlas, atlas_1ds.at( 2 ) );
            bc = std::make_shared<const TPBasisComplex>( atlas, bc, bc_1ds.at( 2 ) );

            return TPSplineSpace( bc, ss, ss_1ds.at( 2 ) );
        }
        else
            throw std::runtime_error( "Unsupported BSpline dimesnion" );
    }

    Eigen::SparseMatrix<double> refinementOp( const TPSplineSpace& coarse, const TPSplineSpace& fine, const double param_tol )
    {
        const auto fine_comps = tensorProductComponentSplines( fine );
        const auto coarse_comps = tensorProductComponentSplines( coarse );

        if( fine_comps.size() != coarse_comps.size() )
            throw std::runtime_error(
                "Cannot create refinement operator between spaces of different parametric dimensions." );

        if( fine_comps.size() == 0 )
            throw std::runtime_error( "Refinement operator only available on tensor products of 1d splines." );

        const util::IndexVec degrees = [&]() {
            util::IndexVec degrees;
            std::transform( fine_comps.begin(),
                            fine_comps.end(),
                            std::back_inserter( degrees ),
                            []( const std::shared_ptr<const BSplineSpace1d>& bspl ) {
                                return bspl->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 );
                            } );

            for( size_t i = 0; i < degrees.size(); i++ )
            {
                if( degrees.at( i ) != coarse_comps.at( i )->basisComplex().defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) )
                    throw std::runtime_error(
                        "p refinement is not supported. Fine and coarse spline spaces must have the same degree." );
            }

            return degrees;
        }();

        SmallVector<KnotVector, 3> kvs_coarse;
        std::transform( coarse_comps.begin(), coarse_comps.end(), std::back_inserter( kvs_coarse ), []( const auto& comp ) { return comp->knotVector(); } );

        SmallVector<KnotVector, 3> kvs_fine;
        std::transform( fine_comps.begin(), fine_comps.end(), std::back_inserter( kvs_fine ), []( const auto& comp ) { return comp->knotVector(); } );

        return refinementOp( kvs_coarse, kvs_fine, degrees, param_tol );
    }
}