#include <BSplineSpace1d.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParametricAtlas.hpp>
#include <BasisComplex1d.hpp>
#include <Eigen/Sparse>
#include <KnotVector.hpp>
#include <BezierTransfer.hpp>

namespace basis
{
    BSplineSpace1d::BSplineSpace1d( const std::shared_ptr<const BasisComplex1d>& bc, const KnotVector& kv )
        : mBasisComplex( bc ), mKnotVector( kv )
    {
        const size_t degree = bc->defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 );
        if( kv.uniqueKnotMultiplicities().at( 0 ).second != degree + 1 or kv.uniqueKnotMultiplicities().back().second != degree + 1 )
            throw std::runtime_error( "Can only create BSplineSpace1d on open knot vector" );

        using SparseMatrixXd = Eigen::SparseMatrix<double>;

        const topology::CombinatorialMap1d& cmap = bc->parametricAtlas().cmap();
        const std::vector<LocalExtraction> local_extractions = localExtractions( kv, degree );

        iterateCellsWhile( cmap, 1, [&]( const topology::Edge& e ) {
            const size_t elem_ii = e.dart().id();
            const LocalExtraction& local_extraction = local_extractions.at( elem_ii );
            mConnectivity.emplace( elem_ii, local_extraction.connectivity );
            mExtractionOps.emplace( elem_ii, local_extraction.extraction );
            return true;
        } );

        /// Add vertex connectivity and operators
        const SparseMatrixXd C = globalExtractionOp( kv, degree );
        // Vertex extraction represents a single shared vertex value and is valid for C >= 0 continuity.
        // For C^{-1} spaces, interior vertices have two one-sided traces that should be handled by the
        // side-aware trace layer rather than this vertex extraction path.
        const auto get_vertex_connectivity = [&degree, &C]( const size_t vertex_ii ) {
            std::vector<FunctionId> rows;
            for( SparseMatrixXd::InnerIterator it( C, vertex_ii * degree ); it; ++it )
                rows.push_back( FunctionId( it.row() ) );
            return rows;
        };

        const auto get_vertex_operator = [&degree, &C]( const size_t vertex_ii, const std::vector<FunctionId>& conn_elem ) {
            Eigen::MatrixXd C_elem = Eigen::MatrixXd::Zero( conn_elem.size(), 1 );
            Eigen::Index row = 0;
            for( SparseMatrixXd::InnerIterator it( C, vertex_ii * degree ); it; ++it )
                C_elem( row++, 0 ) = it.value();
            return C_elem;
        };

        iterateCellsWhile( cmap, 0, [&]( const topology::Vertex& v ) {
            const size_t vertex_ii = v.dart().id();
            const std::vector<FunctionId>& conn_vert = mVertConnectivity.emplace( vertex_ii, get_vertex_connectivity( vertex_ii ) ).first->second;
            mVertExtractionOps.emplace( vertex_ii, get_vertex_operator( vertex_ii, conn_vert ) );
            return true;
        } );

        const size_t end_vert = cellCount( cmap, 0 );
        const std::vector<FunctionId>& conn_vert = mVertConnectivity.emplace( end_vert, get_vertex_connectivity( end_vert ) ).first->second;
        mVertExtractionOps.emplace( end_vert, get_vertex_operator( end_vert, conn_vert ) );
    }

    const BasisComplex1d& BSplineSpace1d::basisComplex() const
    {
        return *mBasisComplex;
    }

    Eigen::MatrixXd BSplineSpace1d::extractionOperator( const topology::Cell& c ) const
    {
        if( c.dim() == 1 )
        {
            return mExtractionOps.at( c.dart().id() );
        }
        else if( c.dim() == 0 )
        {
            return mVertExtractionOps.at( c.dart().id() );
        }
        else
            throw std::runtime_error( "Too large cell dimension for BSplineSpace1d" );
    }

    std::vector<FunctionId> BSplineSpace1d::connectivity( const topology::Cell& c ) const
    {
        if( c.dim() == 1 )
        {
            return mConnectivity.at( c.dart().id() );
        }
        else if( c.dim() == 0 )
        {
            return mVertConnectivity.at( c.dart().id() );
        }
        else
            throw std::runtime_error( "Extraction operators only supported for elements" );
    }

    Eigen::MatrixXd BSplineSpace1d::endVertexExtractionOperator() const
    {
        return mVertExtractionOps.at( cellCount( mBasisComplex->parametricAtlas().cmap(), 0 ) );
    }

    std::vector<FunctionId> BSplineSpace1d::endVertexConnectivity() const
    {
        return mVertConnectivity.at( cellCount( mBasisComplex->parametricAtlas().cmap(), 0 ) );
    }

    size_t BSplineSpace1d::numFunctions() const
    {
        return mKnotVector.size() - mBasisComplex->defaultParentBasis().mBasisGroups.at( 0 ).degrees.at( 0 ) - 1;
    }
}
