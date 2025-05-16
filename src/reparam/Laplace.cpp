#include <Laplace.hpp>
#include <Simplex.hpp>
#include <Logging.hpp>
#include <SimplexUtilities.hpp>
#include <TetMeshCombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CombinatorialMapBoundary.hpp>
#include <Eigen/Sparse>
#include <numeric>
#include <CutCombinatorialMap.hpp>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>

#define LOG_LAPLACE 0

namespace reparam
{
    Timer t( 11 );

    std::vector<double> cotanEdgeWeights3d( const topology::CombinatorialMap& map,
                                            const VertexPositionsFunc& vertex_position,
                                            const std::vector<Normal>& normals )
    {
        const auto edge_ids = indexingOrError( map, 1 );
        const size_t n_edges = cellCount( map, 1 );
        std::vector<double> weights( n_edges, 0 );

        iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
            double& weight = weights.at( edge_ids( e ) );
            iterateAdjacentCells( map, e, 3, [&]( const topology::Volume& v ) {
                const topology::Edge op_edge( phi( map, { 1, 2, -1 }, v.dart() ).value() );
                const double factor = edgeLength( map, vertex_position, op_edge ) / 12;
                weight += factor * dihedralCotangent( map, op_edge, normals );
                return true;
            } );
            return true;
        } );
        return weights;
    }

    std::vector<double> barycentricDualEdgeWeights3d( const topology::CombinatorialMap& map,
                                                      const VertexPositionsFunc& v_positions )
    {
        const auto edge_ids = indexingOrError( map, 1 );
        const size_t n_edges = cellCount( map, 1 );
        std::vector<double> weights( n_edges, 0 );

        iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
            double& weight = weights.at( edge_ids( e ) );

            const Eigen::Vector3d edge_mid = 0.5 * ( v_positions( topology::Vertex( e.dart() ) ) + v_positions( topology::Vertex( phi( map, 1, e.dart() ).value() ) ) );
            const double edge_len = edgeLength( map, v_positions, e );
            const double factor = 6 / ( edge_len * edge_len );
            iterateAdjacentCells( map, e, 3, [&]( const topology::Volume& v ) {
                const Eigen::Vector3d vert_position = v_positions( topology::Vertex( v.dart() ) );
                const Tetrahedron tet = tetOfVolume( map, v_positions, v );
                const Eigen::Vector3d tet_circum = centroid( tet );
                const Eigen::Vector3d face_circum = centroid( triangleOfFace<3>( map, v_positions, topology::Face( v.dart() ) ) );
                const Eigen::Vector3d opp_face_circum = centroid( triangleOfFace<3>( map, v_positions, topology::Face( phi( map, 2, v.dart() ).value() ) ) );
                weight += tetVolume( Tetrahedron( { edge_mid, tet_circum, face_circum, vert_position } ) );
                weight += tetVolume( Tetrahedron( { edge_mid, opp_face_circum, tet_circum, vert_position } ) );

                return true;
            } );

            weight *= factor;

            return true;
        } );

        return weights;
    }

    std::vector<double> voronoiDualEdgeWeights3d( const topology::CombinatorialMap& map, const VertexPositionsFunc& v_positions )
    {
        const auto edge_ids = indexingOrError( map, 1 );
        const size_t n_edges = cellCount( map, 1 );
        std::vector<double> weights( n_edges, 0 );

        iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
            double& weight = weights.at( edge_ids( e ) );

            const Eigen::Vector3d edge_mid = 0.5 * ( v_positions( topology::Vertex( e.dart() ) ) + v_positions( topology::Vertex( phi( map, 1, e.dart() ).value() ) ) );
            const double edge_len = edgeLength( map, v_positions, e );
            const double factor = 6 / ( edge_len * edge_len );
            iterateAdjacentCells( map, e, 3, [&]( const topology::Volume& v ) {
                const Eigen::Vector3d vert_position = v_positions( topology::Vertex( v.dart() ) );
                const Tetrahedron tet = tetOfVolume( map, v_positions, v );
                const Triangle<3> face1 = triangleOfFace<3>( map, v_positions, topology::Face( v.dart() ) );
                const Triangle<3> face2 = triangleOfFace<3>( map, v_positions, topology::Face( phi( map, 2, v.dart() ).value() ) );
                const Eigen::Vector3d tet_circum = circumcenter( tet );
                const Eigen::Vector3d face_circum = circumcenter( face1 );
                const Eigen::Vector3d opp_face_circum = circumcenter( face2 );
                weight += tetVolume( Tetrahedron( { edge_mid, tet_circum, face_circum, vert_position } ) );
                weight += tetVolume( Tetrahedron( { edge_mid, opp_face_circum, tet_circum, vert_position } ) );

                // NOTE: if we check here for circumcenters outside the tetrahedron or the triangles and calculate the
                // barycentric dual weight instead, we get the hybrid weights mentioned in the Alexa et al. paper.

                return true;
            } );

            weight *= factor;
            return true;
        } );
        return weights;
    }

    std::vector<double> inverseLengthEdgeWeights( const topology::CombinatorialMap& map,
                                                  const VertexPositionsFunc& v_positions )
    {
        const auto edge_ids = indexingOrError( map, 1 );
        const size_t n_edges = cellCount( map, 1 );
        std::vector<double> weights( n_edges, 0 );

        iterateCellsWhile( map, 1, [&]( const topology::Edge& e ) {
            weights.at( edge_ids( e ) ) = 1.0 / edgeLength( map, v_positions, e );
            return true;
        } );
        return weights;
    }

    std::vector<double> edgeWeightsLaplace3d( const topology::CombinatorialMap& map,
                                             const VertexPositionsFunc& vertex_position,
                                             const std::vector<Normal>& normals,
                                             const Laplace3dEdgeWeights& edge_weights )
    {
        switch( edge_weights )
        {
            case Laplace3dEdgeWeights::Cotangent:
                return cotanEdgeWeights3d( map, vertex_position, normals );
            case Laplace3dEdgeWeights::VoronoiDual:
                return voronoiDualEdgeWeights3d( map, vertex_position );
            case Laplace3dEdgeWeights::BarycentricDual:
                return barycentricDualEdgeWeights3d( map, vertex_position );
            case Laplace3dEdgeWeights::InverseLength:
                return inverseLengthEdgeWeights( map, vertex_position );
            case Laplace3dEdgeWeights::Uniform:
                return std::vector<double>( cellCount( map, 1 ), 1.0 );
        }
    }

    Eigen::SparseVector<double>
        laplaceOperatorRowSparse( const topology::CombinatorialMap& map,
                                  const topology::Vertex& v1,
                                  const std::function<double( const topology::Edge& )>& edge_weights,
                                  const int n_verts )
    {
        const auto vertex_ids = indexingOrError( map, 0 );
        Eigen::SparseVector<double> out( n_verts );
        out.reserve( 10 ); // FIXME
        const VertexId vid_ref = vertex_ids( v1 );
        t.start( 7 );
        iterateAdjacentCells( map, v1, 1, [&]( const topology::Edge& e ) {
            const double edge_weight = edge_weights( e );
            const VertexId vid1 = vertex_ids( topology::Vertex( e.dart() ) );
            const VertexId vid2 = vertex_ids( topology::Vertex( phi( map, 1, e.dart() ).value() ) );

            if( vid1 == vid_ref )
            {
                out.coeffRef( vid1.id() ) -= edge_weight;
                out.coeffRef( vid2.id() ) += edge_weight;
            }
            else
            {
                out.coeffRef( vid2.id() ) -= edge_weight;
                out.coeffRef( vid1.id() ) += edge_weight;
            }
            return true;
        } );
        t.stop( 7 );

        return out;
    }

    Eigen::VectorXd sweepEmbedding( const topology::TetMeshCombinatorialMap& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals,
                                    const Laplace3dEdgeWeights& edge_weights_type )
    {
        const auto vertex_ids = indexingOrError( map, 0 );
        const auto vertex_position = [&]( const topology::Vertex& v ) {
            return map.simplicialComplex().points.at( vertex_ids( v ) );
        };

        const auto edge_ids = indexingOrError( map, 1 );
        const std::vector<double> edge_weights = edgeWeightsLaplace3d( map, vertex_position, normals, edge_weights_type );
        const auto edge_weights_func = [&]( const topology::Edge& e ) { return edge_weights.at( edge_ids( e ) ); };

        const auto constraints = [&]( const topology::Vertex& v ) -> std::optional<Eigen::VectorXd> {
            if( zero_bcs.at( vertex_ids( v ) ) )
                return Eigen::Matrix<double, 1, 1>( 0.0 );
            else if( one_bcs.at( vertex_ids( v ) ) )
                return Eigen::Matrix<double, 1, 1>( 1.0 );
            else
                return {};
        };

        const size_t n_constraints = std::accumulate( zero_bcs.begin(), zero_bcs.end(), 0 ) +
                                     std::accumulate( one_bcs.begin(), one_bcs.end(), 0 );

        return solveLaplaceSparse( map, edge_weights_func, constraints, n_constraints, 1 );
    }

    Eigen::MatrixX2d
        tutteEmbedding( const topology::CombinatorialMap& map,
                        const VertexPositionsFunc& vert_positions,
                        const std::function<std::optional<Eigen::Vector2d>( const topology::Vertex& )>& constraints,
                        const bool shape_preserving )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Tutte embedding only supports 2d maps" );

        const size_t n_bdry_verts = [&map]() {
            const topology::CombinatorialMapBoundary bdry( map );
            return cellCount( bdry, 0 );
        }();

        const auto constraints_wrapper = [&constraints]( const topology::Vertex& v ) -> std::optional<Eigen::VectorXd> {
            return constraints( v ).transform( [&]( const Eigen::Vector2d& vec ) { return Eigen::VectorXd( vec ); } );
        };

        const auto edge_weights = [&]() -> std::function<double( const topology::Edge& )> {
            if( shape_preserving )
                return [&]( const topology::Edge& e ) { return 1.0 / edgeLength( map, vert_positions, e ); };
            else
                return []( const auto& ) { return 1.0; };
        }();

        return solveLaplaceSparse( map, edge_weights, constraints_wrapper, n_bdry_verts, map.dim() );
    }

    /// See https://dl.acm.org/doi/10.1145/2816795.2818099
    /// "Orbifold Tutte embeddings," by Aigerman and Lipman.  We are using the type (i) orbifold tile.
    Eigen::MatrixX2d tutteOrbifoldEmbedding( const topology::CutCombinatorialMap& map,
                                             const VertexPositionsFunc& vert_positions,
                                             const std::array<topology::Vertex, 3>& cone_vertices,
                                             const bool shape_preserving )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Tutte embedding only supports 2d maps" );

        const auto edge_weights = [&]() -> std::function<double( const topology::Edge& )> {
            if( shape_preserving )
                return [&]( const topology::Edge& e ) { return 1.0 / edgeLength( map, vert_positions, e ); };
            else
                return []( const auto& ) { return 1.0; };
        }();

        const topology::Vertex start_v( maybeBoundaryDart( map, cone_vertices.at( 0 ) ).value() );
        const topology::CombinatorialMapBoundary bdry( map, { start_v.dart() } );
        const size_t n_bdry_verts = cellCount( bdry, 0 );
        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );

        /*
           - For all the interior vertices of map, use the normal laplaceOperatorRowSparse.
           - For the non-cone boundary vertices, use a modified laplaceOperatorRowSparse that
             acts as if both copies of the vertices and their neighborhood are at the same place.
             Also add another equation that constrains them to be the correct rotations of each other.
           - For the cone vertices, use constraints.
        */

        t.start( 0 );

        using SparseVectorXd = Eigen::SparseVector<double>;
        using SparseMatrixXd = Eigen::SparseMatrix<double>;
        std::map<size_t, Eigen::Index> unknown_verts;

        const size_t n_verts = cellCount( map, 0 );

        const auto vertex_ids = indexingOrError( map, 0 );
        const topology::CombinatorialMap& uncut_map = map.baseMap();
        const auto uncut_vertex_ids = indexingOrError( uncut_map, 0 );

        constexpr size_t n_constrained_verts = 4;
        const auto other_side_of_cut = [&]( const topology::Vertex& v ) {
            const topology::Vertex v_bdry = [&](){
                const auto maybe = maybeBoundaryDart( map, v );
                if( not maybe.has_value() )
                {
                    std::cout << "Vertex " << v << " is not on the boundary!" << std::endl;
                    std::cout << vert_positions( v ).transpose() << std::endl;
                }
                return maybe.value();
            }();
            const topology::Vertex other_v_bdry( phi( uncut_map, {2,1}, v_bdry.dart() ).value() );
            return other_v_bdry;
        };
        const std::array<std::pair<size_t, Eigen::Index>, 4> constrained_verts = [&]() -> std::array<std::pair<size_t, Eigen::Index>, 4> {
            const topology::Vertex other_mid_vert = other_side_of_cut( cone_vertices.at( 1 ) );

            return { { { vertex_ids( cone_vertices.at( 0 ) ), 0 },
                       { vertex_ids( cone_vertices.at( 1 ) ), 1 },
                       { vertex_ids( cone_vertices.at( 2 ) ), 2 },
                       { vertex_ids( other_mid_vert ), 3 } } };
        }();

        std::vector<Eigen::Triplet<double>> L_triplets;
        L_triplets.reserve( 4 * cellCount( map, 1 ) + 2 * n_verts + 5 * ( n_bdry_verts - n_constrained_verts ) );

        const Eigen::Matrix<double, 8, 1> BCs = ( Eigen::Matrix<double, 8, 1>() << 0, -1, 1, 0, 0, 1, -1, 0 ).finished();

        const auto add_doubled_row = []( const SparseVectorXd& row, const Eigen::Index i, std::vector<Eigen::Triplet<double>>& triplets ) {
            for( SparseVectorXd::InnerIterator it( row ); it; ++it )
            {
                triplets.emplace_back( 2 * i, 2 * it.index(), it.value() );
                triplets.emplace_back( 2 * i + 1, 2 * it.index() + 1, it.value() );
            }
        };

        const auto add_rotated_doubled_row = [&add_doubled_row]( const SparseVectorXd& row,
                                                                 const SparseMatrixXd& rot,
                                                                 const Eigen::Index i,
                                                                 std::vector<Eigen::Triplet<double>>& triplets_out ) {
            using SparseRowMatrixXd = Eigen::SparseMatrix<double, Eigen::RowMajor>;
            SparseRowMatrixXd doubled_row( 2, 2 * row.size() );
            std::vector<Eigen::Triplet<double>> row_triplets;
            row_triplets.reserve( 2 * row.nonZeros() );
            add_doubled_row( row, 0, row_triplets );
            doubled_row.setFromTriplets( row_triplets.begin(), row_triplets.end() );
            SparseRowMatrixXd rotated_row = rot * doubled_row;
            for( SparseRowMatrixXd::InnerIterator it( rotated_row, 0 ); it; ++it )
                triplets_out.emplace_back( 2 * i, it.index(), it.value() );
            for( SparseRowMatrixXd::InnerIterator it( rotated_row, 1 ); it; ++it )
                triplets_out.emplace_back( 2 * i + 1, it.index(), it.value() );
        };

        const auto add_constrain_rotation_rows = []( const SparseMatrixXd& rot,
                                                     const Eigen::Index row,
                                                     const Eigen::Index i_col,
                                                     const Eigen::Index j_col,
                                                     const Eigen::Index common_vert_col,
                                                     std::vector<Eigen::Triplet<double>>& triplets_out ) {
            triplets_out.emplace_back( row, 2 * i_col, -1 );
            triplets_out.emplace_back( row + 1, 2 * i_col + 1, -1 );

            triplets_out.emplace_back( row, 2 * j_col, rot.coeff( 0, 0 ) );
            triplets_out.emplace_back( row, 2 * j_col + 1, rot.coeff( 0, 1 ) );
            triplets_out.emplace_back( row + 1, 2 * j_col, rot.coeff( 1, 0 ) );
            triplets_out.emplace_back( row + 1, 2 * j_col + 1, rot.coeff( 1, 1 ) );

            triplets_out.emplace_back( row, 2 * common_vert_col, 1 - rot.coeff( 0, 0 ) );
            triplets_out.emplace_back( row, 2 * common_vert_col + 1, -rot.coeff( 0, 1 ) );
            triplets_out.emplace_back( row + 1, 2 * common_vert_col, -rot.coeff( 1, 0 ) );
            triplets_out.emplace_back( row + 1, 2 * common_vert_col + 1, 1 - rot.coeff( 1, 1 ) );
        };

        t.start( 1 );
        iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
            const size_t vid = vertex_ids( v );
            if( vid >= n_verts )
                throw std::runtime_error( "Solving a Laplace system requires a contiguous zero based vertex indexing" );

            if( not boundaryAdjacent( map, v ) )
            {
                t.start( 2 );
                const SparseVectorXd row = laplaceOperatorRowSparse( map, v, edge_weights, n_verts );
                t.stop( 2 );
                const Eigen::Index i = unknown_verts.size();
                add_doubled_row( row, i, L_triplets );
                unknown_verts.emplace( vid, i );
            }
            return true;
        } );
        t.stop( 1 );

        const auto sparse_rotation = []( const double angle ) -> SparseMatrixXd {
            const Eigen::Matrix2d rot = Eigen::Rotation2Dd( angle ).toRotationMatrix();
            SparseMatrixXd sparse_rot( 2, 2 );
            if( std::abs( rot( 0, 0 ) ) > 1e-12 ) sparse_rot.insert( 0, 0 ) = rot( 0, 0 );
            if( std::abs( rot( 0, 1 ) ) > 1e-12 ) sparse_rot.insert( 0, 1 ) = rot( 0, 1 );
            if( std::abs( rot( 1, 0 ) ) > 1e-12 ) sparse_rot.insert( 1, 0 ) = rot( 1, 0 );
            if( std::abs( rot( 1, 1 ) ) > 1e-12 ) sparse_rot.insert( 1, 1 ) = rot( 1, 1 );
            sparse_rot.makeCompressed();
            return sparse_rot;
        };

        t.start( 10 );
        {
            // We take phi1s, add the rows for the vertices on both sides of the boundary, repeat until
            // we get to the next cone vertex, skip that, and continue to the last cone vertex.
            topology::Dart curr_d = start_v.dart();
            SparseMatrixXd rot = sparse_rotation( -0.5 * std::numbers::pi );
            while( bdry_vertex_ids( topology::Vertex( curr_d ) ) != constrained_verts.at( 1 ).first )
            {
                const topology::Vertex curr_v( bdry.toUnderlyingCell( topology::Vertex( curr_d ) ) );
                const topology::Vertex other_v = other_side_of_cut( curr_v );
                const Eigen::Index i = unknown_verts.size();

                t.start( 2 );
                const SparseVectorXd row_1 = laplaceOperatorRowSparse( map, curr_v, edge_weights, n_verts );
                const SparseVectorXd row_2 = laplaceOperatorRowSparse( map, other_v, edge_weights, n_verts );
                t.stop( 2 );

                add_doubled_row( row_1, i, L_triplets );
                add_rotated_doubled_row( row_2, rot, i, L_triplets );
                add_constrain_rotation_rows( rot, 2 * i + 2, vertex_ids( curr_v ), vertex_ids( other_v ), constrained_verts.at( 0 ).first, L_triplets );

                unknown_verts.emplace( vertex_ids( curr_v ), i );
                unknown_verts.emplace( vertex_ids( other_v ), i + 1 );
                curr_d = phi( bdry, -1, curr_d ).value();
            }

            rot = sparse_rotation( 0.5 * std::numbers::pi );
            curr_d = phi( bdry, -1, curr_d ).value();

            while( bdry_vertex_ids( topology::Vertex( curr_d ) ) != constrained_verts.at( 2 ).first )
            {
                const topology::Vertex curr_v( bdry.toUnderlyingCell( topology::Vertex( curr_d ) ) );
                const topology::Vertex other_v = other_side_of_cut( curr_v );
                const Eigen::Index i = unknown_verts.size();

                t.start( 2 );
                const SparseVectorXd row_1 = laplaceOperatorRowSparse( map, curr_v, edge_weights, n_verts );
                const SparseVectorXd row_2 = laplaceOperatorRowSparse( map, other_v, edge_weights, n_verts );
                t.stop( 2 );

                add_doubled_row( row_1, i, L_triplets );
                add_rotated_doubled_row( row_2, rot, i, L_triplets );
                add_constrain_rotation_rows( rot, 2 * i + 2, vertex_ids( curr_v ), vertex_ids( other_v ), constrained_verts.at( 2 ).first, L_triplets );

                unknown_verts.emplace( vertex_ids( curr_v ), i );
                unknown_verts.emplace( vertex_ids( other_v ), i + 1 );

                curr_d = phi( bdry, -1, curr_d ).value();
            }
        }
        t.stop( 10 );

        t.start( 3 );
        std::vector<Eigen::Triplet<double>> L_II_triplets;
        L_II_triplets.reserve( L_triplets.size() );
        std::vector<Eigen::Triplet<double>> L_IB_triplets;
        L_IB_triplets.reserve( L_triplets.size() );
        for( const auto& t : L_triplets )
        {
            const size_t vert_id = t.col() / 2;
            const auto find_it = std::find_if( constrained_verts.begin(), constrained_verts.end(), [&]( const auto& pr ) { return pr.first == (size_t)vert_id; } );
            if( find_it != constrained_verts.end() )
            {
                L_IB_triplets.emplace_back( t.row(), 2 * find_it->second + t.col() % 2, t.value() );
            }
            else
            {
                L_II_triplets.emplace_back( t.row(), 2 * unknown_verts.at( vert_id ) + t.col() % 2, t.value() );
            }
        }
        SparseMatrixXd L_II( 2 * ( n_verts - n_constrained_verts), 2 * ( n_verts - n_constrained_verts ) );
        L_II.setFromTriplets( L_II_triplets.begin(), L_II_triplets.end() );

        SparseMatrixXd L_IB( 2 * ( n_verts - n_constrained_verts ), 2 * n_constrained_verts );
        L_IB.setFromTriplets( L_IB_triplets.begin(), L_IB_triplets.end() );

        LOG( LOG_LAPLACE ) << "BCs: " << std::endl << BCs << std::endl << std::endl;

        const Eigen::MatrixXd rhs = -L_IB * BCs;
        t.stop( 3 );

        LOG( LOG_LAPLACE ) << "L_II:\n" << Eigen::MatrixXd( L_II ) << std::endl << std::endl;
        LOG( LOG_LAPLACE ) << "L_IB:\n" << Eigen::MatrixXd( L_IB ) << std::endl << std::endl;
        LOG( LOG_LAPLACE ) << "rhs:\n" << rhs.transpose() << std::endl << std::endl;

        LOG( LOG_LAPLACE ) << "About to solve\n";

        t.start( 4 );
        const Eigen::MatrixXd ans = [&]() -> Eigen::MatrixXd {
            Eigen::SparseLU<SparseMatrixXd> solver( L_II );
            return solver.solve( rhs );
        }();

        t.stop( 4 );

        LOG( LOG_LAPLACE ) << "Solved\n" << ans.transpose() << std::endl;

        LOG( LOG_LAPLACE ) << "Assembling result\n";

        t.start( 5 );
        Eigen::MatrixX2d result( n_verts, 2 );
        for( const auto& pr : unknown_verts )
        {
            result( pr.first, 0 ) = ans( 2 * pr.second );
            result( pr.first, 1 ) = ans( 2 * pr.second + 1 );
        }
        for( const auto& pr : constrained_verts )
        {
            result( pr.first, 0 ) = BCs( 2 * pr.second );
            result( pr.first, 1 ) = BCs( 2 * pr.second + 1 );
        }
        t.stop( 5 );
        t.stop( 0 );

        std::cout << result << std::endl;

        LOG( LOG_LAPLACE ) << "returning result\n";

        LOG( LOG_LAPLACE ) << "Total time: " << t.stop( 0 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Weights time: " << t.stop( 9 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| | Edge length time: " << t.stop( 8 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| | cot time: " << t.stop( 6 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Loop time: " << t.stop( 1 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Boundary time: " << t.stop( 10 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| | Row time: " << t.stop( 2 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| | | Loop time: " << t.stop( 7 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Assembly time: " << t.stop( 3 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Solve time: " << t.stop( 4 ) << std::endl;
        LOG( LOG_LAPLACE ) << "| Format time: " << t.stop( 5 ) << std::endl;

        return result;
        return Eigen::MatrixX2d();
    }

    Eigen::MatrixXd
        solveLaplaceSparse( const topology::CombinatorialMap& map,
                            const std::function<double( const topology::Edge& )>& edge_weights,
                            const std::function<std::optional<Eigen::VectorXd>( const topology::Vertex& )>& constraints,
                            const size_t n_constrained_verts,
                            const size_t data_dim )
    {
        t.start( 0 );

        using SparseVectorXd = Eigen::SparseVector<double>;
        using SparseMatrixXd = Eigen::SparseMatrix<double>;
        std::map<size_t, Eigen::Index> interior_verts;
        std::map<size_t, Eigen::Index> boundary_verts;

        const size_t n_verts = cellCount( map, 0 );

        const auto vertex_ids = indexingOrError( map, 0 );

        std::vector<Eigen::Triplet<double>> L_triplets;
        L_triplets.reserve( 2 * cellCount( map, 1 ) + n_verts );

        Eigen::MatrixXd BCs( n_constrained_verts, data_dim );

        t.start( 1 );
        iterateCellsWhile( map, 0, [&]( const topology::Vertex& v ) {
            const size_t vid = vertex_ids( v );
            if( vid >= n_verts )
                throw std::runtime_error( "Solving a Laplace system requires a contiguous zero based vertex indexing" );
            const auto maybe_constrained = constraints( v );
            if( maybe_constrained.has_value() )
            {
                const Eigen::Index i = boundary_verts.size();
                BCs.row( i ) = maybe_constrained.value();
                boundary_verts.emplace( vid, i );
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
                interior_verts.emplace( vid, i );
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
        SparseMatrixXd L_II( n_verts - n_constrained_verts, n_verts - n_constrained_verts );
        L_II.setFromTriplets( L_II_triplets.begin(), L_II_triplets.end() );

        std::vector<Eigen::Triplet<double>> L_IB_triplets;
        L_IB_triplets.reserve( L_triplets.size() );
        for( const auto& t : L_triplets )
        {
            const auto find_it = boundary_verts.find( t.col() );
            if( find_it != boundary_verts.end() ) L_IB_triplets.emplace_back( t.row(), find_it->second, t.value() );
        }
        SparseMatrixXd L_IB( n_verts - n_constrained_verts, n_constrained_verts );
        L_IB.setFromTriplets( L_IB_triplets.begin(), L_IB_triplets.end() );

        LOG( LOG_LAPLACE ) << "BCs: " << std::endl << BCs << std::endl << std::endl;

        const Eigen::MatrixXd rhs = -L_IB * BCs;
        t.stop( 3 );

        LOG( LOG_LAPLACE ) << "L_II:\n" << Eigen::MatrixXd( L_II ) << std::endl << std::endl;
        LOG( LOG_LAPLACE ) << "L_IB:\n" << Eigen::MatrixXd( L_IB ) << std::endl << std::endl;
        LOG( LOG_LAPLACE ) << "rhs:\n" << rhs.transpose() << std::endl << std::endl;

        LOG( LOG_LAPLACE ) << "About to solve\n";

        t.start( 4 );
        const Eigen::MatrixXd ans = [&]() -> Eigen::MatrixXd {
            if( map.dim() == 3 )
            {
                Eigen::ConjugateGradient<SparseMatrixXd, Eigen::Lower | Eigen::Upper> solver( L_II );
                return solver.solve( rhs );
            }
            else
            {
                Eigen::SimplicialLDLT<SparseMatrixXd> solver( L_II );
                return solver.solve( rhs );
            }
        }();

        t.stop( 4 );

        LOG( LOG_LAPLACE ) << "Assembling result\n";

        t.start( 5 );
        Eigen::MatrixXd result( n_verts, data_dim );
        for( const auto& pr : interior_verts ) result.row( pr.first ) = ans.row( pr.second );
        for( const auto& pr : boundary_verts ) result.row( pr.first ) = BCs.row( pr.second );
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
} // namespace reparam