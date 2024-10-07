#include <LeastSquaresFitting.hpp>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

namespace fitting
{
    Eigen::MatrixXd leastSquaresFitting(
        eval::SplineSpaceEvaluator& evaler,
        const size_t n_points,
        const size_t spatial_dim,
        const std::function<void( const std::function<void(const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd&)>& )>& point_iterator )
    {
        std::vector<Eigen::Triplet<double>> triplets;//FIXME: reserve

        Eigen::MatrixXd rhs( n_points, spatial_dim );

        size_t row = 0;
        point_iterator( [&]( const topology::Cell& elem, const param::ParentPoint& ppt, const Eigen::Vector3d& field_pt ) {
            evaler.localizeElement( elem ); // TODO: Move this and the connectivity out into a less hot part of the loop
            evaler.localizePoint( ppt );
            const std::vector<basis::FunctionId> conn = evaler.splineSpace().connectivity( elem );

            const Eigen::VectorXd basis = evaler.evaluateBasis();

            for( size_t func_ii = 0; func_ii < conn.size(); func_ii++ )
            {
                triplets.emplace_back( row, conn.at( func_ii ), basis( func_ii ) );
            }

            rhs.row( row++ ) = field_pt.transpose();
        } );

        Eigen::SparseMatrix<double> A( n_points, evaler.splineSpace().numFunctions() );
        A.setFromTriplets( triplets.begin(), triplets.end() );

        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute( A );
        return solver.solve( rhs );
    }
}