#include <Laplace.hpp>
#include <Simplex.hpp>
#include <Logging.hpp>
#include <SimplexUtilities.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <Eigen/Sparse>

#define LOG_LAPLACE 0

Timer t;

double edgeWeight( const topology::TetMeshCombinatorialMap& map, const topology::Edge& e, const std::vector<Normal>& normals )
{
    double weight = 0;
    t.start( 8 );

    t.stop( 8 );
    iterateAdjacentCells( map, e, 3, [&]( const topology::Volume& v ) {
        t.start( 6 );
        const topology::Edge op_edge( phi( map, { 1, 2, -1 }, v.dart() ).value() );
        const double factor = edgeLength( map, op_edge ) / 12;
        weight += factor * dihedralCotangent( map, op_edge, normals );
        t.stop( 6 );
        return true;
    } );

    return weight;
}

std::vector<double> edgeWeights( const topology::TetMeshCombinatorialMap& map, const std::vector<Normal>& normals )
{
    const size_t n_edges = cellCount( map, 1 );
    std::vector<double> weights( n_edges, 0 );

    iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
        weights.at( map.edgeId( e ) ) = edgeWeight( map, e, normals );
        return true;
    } );
    return weights;
}

Eigen::SparseVector<double> laplaceOperatorRowSparse( const topology::TetMeshCombinatorialMap& map,
                                                      const topology::Vertex& v1,
                                                      const std::vector<double>& edge_weights,
                                                      const int n_verts )
{
    Eigen::SparseVector<double> out( n_verts );
    out.reserve( 10 ); // FIXME
    const VertexId vid1 = map.vertexId( v1 );
    t.start( 7 );
    iterateAdjacentCells( map, v1, 1, [&]( const topology::Edge& e ) {
        const double edge_weight = edge_weights.at( map.edgeId( e ) );
        const VertexId vid2 = map.vertexId( topology::Vertex( phi( map, 1, e.dart() ).value() ) );

        out.coeffRef( vid1.id() ) -= edge_weight;
        out.coeffRef( vid2.id() ) += edge_weight;
        return true;
    } );
    t.stop( 7 );

    return out;
}

Eigen::VectorXd solveLaplaceSparse( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals )
{
    t.start( 0 );

    using SparseVectorXd = Eigen::SparseVector<double>;
    using SparseMatrixXd = Eigen::SparseMatrix<double>;
    std::map<Eigen::Index, Eigen::Index> interior_verts;
    std::map<Eigen::Index, Eigen::Index> boundary_verts;

    const size_t n_verts = cellCount( map, 0 );
    const size_t n_bcs = std::accumulate( zero_bcs.begin(), zero_bcs.end(), 0 ) +
                         std::accumulate( one_bcs.begin(), one_bcs.end(), 0 );

    std::vector<Eigen::Triplet<double>> L_triplets;
    L_triplets.reserve( 2 * cellCount( map, 1 ) + n_verts );

    SparseVectorXd BCs( n_bcs );
    BCs.reserve( one_bcs.size() );

    LOG( LOG_LAPLACE ) << "zeros: " << zero_bcs << std::endl;
    LOG( LOG_LAPLACE ) << "ones: " << one_bcs << std::endl;

    t.start( 9 );
    const std::vector<double> edge_weights = edgeWeights( map, normals );
    t.stop( 9 );

    t.start( 1 );
    iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
        const VertexId vid = map.vertexId( v );
        if( zero_bcs.at( vid.id() ) )
        {
            const Eigen::Index i = boundary_verts.size();
            boundary_verts.emplace( vid.id(), i );
        }
        else if( one_bcs.at( vid.id() ) )
        {
            const Eigen::Index i = boundary_verts.size();
            BCs.insert( i ) = 1.0;
            boundary_verts.emplace( vid.id(), i );
        }
        else
        {
            t.start( 2 );
            const SparseVectorXd row = laplaceOperatorRowSparse( map, v, edge_weights, n_verts );
            t.stop( 2 );
            const Eigen::Index i = interior_verts.size();
            for( SparseVectorXd::InnerIterator it( row ); it; ++it )
            {
                L_triplets.emplace_back( i, it.row(), it.value() );
            }
            interior_verts.emplace( vid.id(), i );
        }
        return true;
    } );
    t.stop( 1 );

    t.start( 3 );
    std::vector<Eigen::Triplet<double>> L_II_triplets;
    L_II_triplets.reserve( L_triplets.size() );
    for( const auto& t : L_triplets )
    {
        const auto find_it = interior_verts.find( t.col() );
        if( find_it != interior_verts.end() )
        {
            L_II_triplets.emplace_back( t.row(), find_it->second, t.value() );
        }
    }
    SparseMatrixXd L_II( n_verts - n_bcs, n_verts - n_bcs );
    L_II.setFromTriplets( L_II_triplets.begin(), L_II_triplets.end() );

    std::vector<Eigen::Triplet<double>> L_IB_triplets;
    L_IB_triplets.reserve( L_triplets.size() );
    for( const auto& t : L_triplets )
    {
        const auto find_it = boundary_verts.find( t.col() );
        if( find_it != boundary_verts.end() ) L_IB_triplets.emplace_back( t.row(), find_it->second, t.value() );
    }
    SparseMatrixXd L_IB( n_verts - n_bcs, n_bcs );
    L_IB.setFromTriplets( L_IB_triplets.begin(), L_IB_triplets.end() );

    LOG( LOG_LAPLACE ) << "BCs: " << std::endl << BCs << std::endl << std::endl;

    const SparseVectorXd rhs = -L_IB * BCs;
    t.stop( 3 );

    LOG( LOG_LAPLACE ) << "L_II:\n" << Eigen::MatrixXd( L_II ) << std::endl << std::endl;
    LOG( LOG_LAPLACE ) << "rhs:\n" << Eigen::VectorXd( rhs ).transpose() << std::endl << std::endl;

    LOG( LOG_LAPLACE ) << "About to solve\n";

    t.start( 4 );
    Eigen::ConjugateGradient<SparseMatrixXd, Eigen::Lower | Eigen::Upper> solver( L_II );
    const Eigen::VectorXd ans = solver.solve( rhs );
    t.stop( 4 );

    LOG( LOG_LAPLACE ) << "Assembling result\n";

    t.start( 5 );
    Eigen::VectorXd result( n_verts );
    for( const auto& pr : interior_verts ) result( pr.first ) = ans.coeffRef( pr.second );
    for( const auto& pr : boundary_verts ) result( pr.first ) = BCs.coeffRef( pr.second );
    t.stop( 5 );
    t.stop( 0 );

    LOG( LOG_LAPLACE ) << "returning result\n";

    LOG( LOG_LAPLACE ) << "Total time: " << t.stop( 0 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| Weights time: " << t.stop( 9 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| | Edge length time: " << t.stop( 8 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| | cot time: " << t.stop( 6 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| Loop time: " << t.stop( 1 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| | Row time: " << t.stop( 2 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| | | Loop time: " << t.stop( 7 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| Assembly time: " << t.stop( 3 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| Solve time: " << t.stop( 4 ) << std::endl;
    LOG( LOG_LAPLACE ) << "| Format time: " << t.stop( 5 ) << std::endl;

    return result;
}
