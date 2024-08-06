#include <TPSplineSpace.hpp>
#include <unsupported/Eigen/KroneckerProduct>

namespace basis
{
    TPSplineSpace::TPSplineSpace( const TPBasisComplex& bc, const SplineSpace& source, const BSplineSpace1d& line ) :
        mBasisComplex( bc ), mSource( source ), mLine( line )
    {
        if( &source.basisComplex() != &mBasisComplex.sourceComplex() or &line.basisComplex() != &mBasisComplex.lineComplex() )
            throw std::runtime_error( "Invalid inputs to TPSplineSpace::TPSplineSpace()" );
    }

    const TPBasisComplex& TPSplineSpace::basisComplex() const
    {
        return mBasisComplex;
    }

    Eigen::MatrixXd TPSplineSpace::extractionOperator( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() ) throw std::runtime_error( "Bad cell dimension for spline space" );

        const auto dart_pr = mBasisComplex.parametricAtlas().cmap().unflatten( c.dart() );
        const topology::Cell source_c( std::get<0>( dart_pr ), mSource.basisComplex().parametricAtlas().cmap().dim() );
        const topology::Cell line_c( std::get<1>( dart_pr ), mLine.basisComplex().parametricAtlas().cmap().dim() );
        return Eigen::kroneckerProduct( mLine.extractionOperator( line_c ), mSource.extractionOperator( source_c ) );
    }

    FunctionId TPSplineSpace::flatten( const FunctionId& source_fid, const FunctionId& line_fid ) const
    {
        return FunctionId( source_fid + mSource.numFunctions() * line_fid );
    }


    std::vector<FunctionId> TPSplineSpace::connectivity( const topology::Cell& c ) const
    {
        if( c.dim() != mBasisComplex.parametricAtlas().cmap().dim() ) throw std::runtime_error( "Bad cell dimension for spline space" );

        const auto dart_pr = mBasisComplex.parametricAtlas().cmap().unflatten( c.dart() );
        const std::vector<FunctionId> source_conn = mSource.connectivity(
            topology::Cell( std::get<0>( dart_pr ), mSource.basisComplex().parametricAtlas().cmap().dim() ) );
        const std::vector<FunctionId> line_conn = mLine.connectivity(
            topology::Cell( std::get<1>( dart_pr ), mLine.basisComplex().parametricAtlas().cmap().dim() ) );

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
        return mSource.numFunctions() * mLine.numFunctions();
    }

    std::optional<std::vector<std::reference_wrapper<const BSplineSpace1d>>> constituentSplines( const TPSplineSpace& ss )
    {
        const size_t dim = ss.basisComplex().parametricAtlas().cmap().dim();
        std::vector<std::reference_wrapper<const BSplineSpace1d>> out;
        out.reserve( dim );

        try
        {
            if( dim == 3 )
            {
                const TPSplineSpace& source_primal = dynamic_cast<const TPSplineSpace&>( ss.source() );
                out.push_back( std::cref( dynamic_cast<const BSplineSpace1d&>( source_primal.source() ) ) );
                out.push_back( std::cref( source_primal.line() ) );
            }
            else
                out.push_back( std::cref( dynamic_cast<const BSplineSpace1d&>( ss.source() ) ) );
        }
        catch(const std::bad_cast& e)
        {
            return std::nullopt;
        }
        out.push_back( std::cref( ss.line() ) );
        return out;
    }
}