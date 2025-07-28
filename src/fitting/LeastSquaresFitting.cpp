#include <LeastSquaresFitting.hpp>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <MultiPatchSplineSpace.hpp>
#include <TPSplineSpace.hpp>
#include <set>
#include <SimplicialComplex.hpp>

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

        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.compute( A );
        return solver.solve( rhs );
    }

    Eigen::MatrixXd fitToManifold( const basis::MultiPatchSplineSpace& manifold_ss,
                                   const Eigen::MatrixXd& manifold_cpts,
                                   const basis::MultiPatchSplineSpace& new_ss )
    {
        if( &manifold_ss.basisComplex().parametricAtlas() != &new_ss.basisComplex().parametricAtlas() )
        {
            throw std::runtime_error( "Manifold and new spline space must have the same parametric domain." );
        }

        const Eigen::MatrixXd manifold_cpts_transposed = manifold_cpts.transpose(); // FIXME

        // Get the greville points in the parent domain
        std::vector<std::vector<std::pair<topology::Cell, param::ParentPoint>>> greville_ppts;
        greville_ppts.reserve( new_ss.subSpaces().size() );

        for( size_t patch_ii = 0; patch_ii < new_ss.subSpaces().size(); patch_ii++ )
        {
            const auto& subspace = new_ss.subSpaces().at( patch_ii );
            const basis::TPSplineSpace& tp_subspace = static_cast<const basis::TPSplineSpace&>( *subspace );
            const topology::TPCombinatorialMap& source_cmap = static_cast<const topology::TPCombinatorialMap&>(
                tp_subspace.source().basisComplex().parametricAtlas().cmap() );
            const auto ss_1ds = tensorProductComponentSplines( tp_subspace );

            auto& patch_greville_ppts = greville_ppts.emplace_back();

            SmallVector<std::vector<std::pair<topology::Cell, param::ParentPoint>>, 3> ppts;
            for( const auto& component : ss_1ds )
            {
                const Eigen::VectorXd greville_pts = grevillePoints(
                    component->knotVector(), degrees( component->basisComplex().defaultParentBasis() ).at( 0 ) );
                const std::vector<double> greville_vec( greville_pts.data(), greville_pts.data() + greville_pts.size() );
                ppts.push_back(
                    parentPointsOfParamPoints( greville_vec, component->basisComplex().parametricAtlas(), 1e-9 ) );
            }

            for( const auto& pt3 : ppts.at( 2 ) )
            {
                for( const auto& pt2 : ppts.at( 1 ) )
                {
                    for( const auto& pt1 : ppts.at( 0 ) )
                    {
                        const topology::Volume patch_cell( tp_subspace.basisComplex().parametricAtlas().cmap().flatten(
                            source_cmap.flatten(
                                pt1.first.dart(), pt2.first.dart(), topology::TPCombinatorialMap::TPDartPos::DartPos0 ),
                            pt3.first.dart(),
                            topology::TPCombinatorialMap::TPDartPos::DartPos0 ) );
                        const topology::Volume cell( new_ss.basisComplex().parametricAtlas().cmap().toGlobalDart(
                            patch_ii, patch_cell.dart() ) );
                        patch_greville_ppts.emplace_back( std::make_pair(
                            cell,
                            param::tensorProduct( param::tensorProduct( pt1.second, pt2.second ), pt3.second ) ) );
                    }
                }
            }
        }

        eval::SplineSpaceEvaluator new_evaler( new_ss, 0 );
        eval::SplineSpaceEvaluator manifold_evaler( manifold_ss, 0 );
        return leastSquaresFitting(
            new_evaler,
            new_ss.numFunctions(),
            manifold_cpts.cols(),
            [&]( const std::function<void( const topology::Cell&, const param::ParentPoint&, const Eigen::VectorXd& )>&
                     add_least_squares_point ) {
                const auto& mp_fid_map = new_ss.functionIdMap();
                std::set<basis::FunctionId> seen_fids;
                for( size_t patch_ii = 0; patch_ii < mp_fid_map.size(); patch_ii++ )
                {
                    for( size_t patch_fid = 0; patch_fid < mp_fid_map.at( patch_ii ).size(); patch_fid++ )
                    {
                        const basis::FunctionId fid = mp_fid_map.at( patch_ii ).at( patch_fid );
                        if( seen_fids.count( fid ) > 0 ) continue;
                        seen_fids.insert( fid );

                        const auto& cell_pt = greville_ppts.at( patch_ii ).at( patch_fid );

                        const topology::Cell& cell = cell_pt.first;
                        const param::ParentPoint& ppt = cell_pt.second;

                        manifold_evaler.localizeElement( cell );
                        manifold_evaler.localizePoint( ppt );

                        const Eigen::VectorXd manifold_pt =
                            manifold_evaler.evaluateManifold( manifold_cpts_transposed );

                        add_least_squares_point( cell, ppt, manifold_pt );
                    }
                }
            } );
    }
}