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
#include <CommonUtils.hpp>

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

    double cotanEdgeWeights2d( const topology::CombinatorialMap& map, const VertexPositionsFunc& vertex_position, const topology::Edge& e )
    {
        const Triangle<3> face = triangleOfFace<3>( map, vertex_position, topology::Face( e.dart() ) );
        const auto maybe_phi = phi( map, {2, -1}, e.dart() );
        const Eigen::Vector3d e1 = face.v1 - face.v3;
        const Eigen::Vector3d e2 = face.v2 - face.v3;
        const double cos1 = e1.dot( e2 ) / ( e1.norm() * e2.norm() );
        if( not maybe_phi.has_value() )
        {
            return cos1 / std::sqrt( 1 - cos1 * cos1 );
        }
        const Eigen::Vector3d v4 = vertex_position( topology::Vertex( phi( map, {2, -1}, e.dart() ).value() ) );
        const Eigen::Vector3d e3 = face.v1 - v4;
        const Eigen::Vector3d e4 = face.v2 - v4;
        const double cos2 = e4.dot( e3 ) / ( e4.norm() * e3.norm() );

        return cos1 / std::sqrt( 1 - cos1 * cos1 ) + cos2 / std::sqrt( 1 - cos2 * cos2 );
    }

    double barycentricDualWeights2d( const topology::CombinatorialMap& map, const VertexPositionsFunc& vertex_position, const topology::Edge& e )
    {
        const Triangle<3> face = triangleOfFace<3>( map, vertex_position, topology::Face( e.dart() ) );
        const Eigen::Vector3d v4 = vertex_position( topology::Vertex( phi( map, {2, -1}, e.dart() ).value() ) );
        const Eigen::Vector3d e_mid = 0.5 * ( face.v1 + face.v2 );
        return ( ( face.v3 - e_mid ).norm() + ( v4 - e_mid ).norm() ) / 3.0;
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

    double edgeWeightLaplace2d( const topology::CombinatorialMap& map,
                                const VertexPositionsFunc& vertex_position,
                                const Laplace2dEdgeWeights& edge_weights,
                                const topology::Edge& e )
    {
        switch( edge_weights )
        {
            case Laplace2dEdgeWeights::Cotangent: return cotanEdgeWeights2d( map, vertex_position, e );
            case Laplace2dEdgeWeights::InverseLength: return 1.0 / edgeLength( map, vertex_position, e );
            case Laplace2dEdgeWeights::BarycentricDual: return barycentricDualWeights2d( map, vertex_position, e );
            case Laplace2dEdgeWeights::Uniform: return 1.0;
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
                        const Laplace2dEdgeWeights& edge_weights_type )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Tutte embedding only supports 2d maps" );

        const size_t n_bdry_verts = [&map]() {
            const topology::CombinatorialMapBoundary bdry( map );
            return cellCount( bdry, 0 );
        }();

        const auto constraints_wrapper = [&constraints]( const topology::Vertex& v ) -> std::optional<Eigen::VectorXd> {
            return constraints( v ).transform( [&]( const Eigen::Vector2d& vec ) { return Eigen::VectorXd( vec ); } );
        };

        const auto edge_weights = [&]( const topology::Edge& e ) -> double {
            return edgeWeightLaplace2d( map, vert_positions, edge_weights_type, e );
        };

        return solveLaplaceSparse( map, edge_weights, constraints_wrapper, n_bdry_verts, map.dim() );
    }

    /// See https://dl.acm.org/doi/10.1145/2816795.2818099
    /// "Orbifold Tutte embeddings," by Aigerman and Lipman.  We are using the type (i) orbifold tile.
    Eigen::MatrixX2d tutteOrbifoldEmbedding( const topology::CutCombinatorialMap& map,
                                             const VertexPositionsFunc& vert_positions,
                                             const std::array<topology::Vertex, 3>& cone_vertices,
                                             const Laplace2dEdgeWeights& edge_weights_type )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Tutte embedding only supports 2d maps" );

        const auto edge_weights = [&]( const topology::Edge& e ) -> double {
            return edgeWeightLaplace2d( map, vert_positions, edge_weights_type, e );
        };

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
        std::array<std::pair<size_t, Eigen::Index>, 4> constrained_verts = [&]() -> std::array<std::pair<size_t, Eigen::Index>, 4> {
            const topology::Vertex other_mid_vert = other_side_of_cut( cone_vertices.at( 1 ) );

            return { { { vertex_ids( cone_vertices.at( 0 ) ), 0 },
                       { vertex_ids( cone_vertices.at( 1 ) ), 1 },
                       { vertex_ids( cone_vertices.at( 2 ) ), 2 },
                       { vertex_ids( other_mid_vert ), 3 } } };
        }();

        std::vector<Eigen::Triplet<double>> L_triplets;
        L_triplets.reserve( 4 * cellCount( map, 1 ) + 2 * n_verts + 5 * ( n_bdry_verts - n_constrained_verts ) );

        const Eigen::Matrix<double, 8, 1> BCs = ( Eigen::Matrix<double, 8, 1>() << 0, 0, 1, 0, 1, 1, 0, 1 ).finished();

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
                if( bdry_vertex_ids( topology::Vertex( curr_d ) ) == constrained_verts.at( 3 ).first )
                {
                    std::swap( constrained_verts.at( 1 ).first, constrained_verts.at( 3 ).first );
                    break;
                }
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

    std::pair<Eigen::VectorXd, std::vector<double>> cutTutteBCsAndAngles( const topology::CutCombinatorialMap& cut_cmap,
                                 const VertexPositionsFunc& vert_positions,
                                 const std::map<size_t, Eigen::Index>& constrained_verts,
                                 const std::vector<std::pair<topology::Vertex, topology::Vertex>>& cut_extremities )
    {
        Eigen::VectorXd out( 2 * constrained_verts.size() );
        std::vector<double> angles;
        angles.reserve( cut_extremities.size() );
        const auto& cmap = cut_cmap.baseMap();
        const auto vert_ids = indexingOrError( cmap, 0 );
        const auto cut_vert_ids = indexingOrError( cut_cmap, 0 );

        const auto is_cut_extremity = [&]( const topology::Vertex& v ) -> bool {
            return std::any_of( cut_extremities.begin(), cut_extremities.end(), [&]( const auto& cut_ends ) {
                return onSameVertex( cmap, v.dart(), cut_ends.first.dart() ) or onSameVertex( cmap, v.dart(), cut_ends.second.dart() );
            } );
        };

        const auto is_start_vert = [&, vert_id = vert_ids( cut_extremities.front().first ) ]( const topology::Vertex& v ) {
            return vert_ids( v ) == vert_id and onBoundary( cmap, phi( cmap, -1, v.dart() ).value() );
        };
        const std::map<size_t, Eigen::Vector2d> bdry_constraints =
            reparam::boundaryConstraints( cut_cmap, vert_positions, cut_extremities.size(), is_cut_extremity, is_start_vert );

        const size_t n_sides = cut_extremities.size() * 4;

        for( const auto& pr : bdry_constraints )
        {
            const auto find_it = constrained_verts.find( pr.first );
            if( find_it != constrained_verts.end() )
            {
                const Eigen::Index i = find_it->second;;
                out( 2 * i ) = pr.second( 0 );
                out( 2 * i + 1 ) = pr.second( 1 );
            }
        }

        const auto furthest_reachable = []( const topology::CombinatorialMap& map, const std::vector<int>& ops, const topology::Dart start_d ) {
            topology::Dart current_d = start_d;
            topology::Dart next_d;
            do
            {
                next_d = current_d;
                const auto maybe_phi = phi( map, ops, next_d );
                if( not maybe_phi.has_value() ) return topology::Vertex( current_d );
                next_d = maybe_phi.value();
                current_d = next_d;
            } while( true );
            return topology::Vertex( current_d );
        };

        for( const auto& cut_pair : cut_extremities )
        {
            const Eigen::Vector2d pos1 = bdry_constraints.at( cut_vert_ids( cut_pair.first ) );
            const Eigen::Vector2d pos2 = bdry_constraints.at( cut_vert_ids( cut_pair.second ) );

            const Eigen::Vector2d pos1_ = bdry_constraints.at( cut_vert_ids( furthest_reachable( cmap, {2,1}, cut_pair.first.dart() ) ) );
            const Eigen::Vector2d pos2_ = bdry_constraints.at( cut_vert_ids( furthest_reachable( cmap, {-1,2}, cut_pair.second.dart() ) ) );

            angles.push_back( - std::atan2( pos2_( 1 ) - pos1_( 1 ), pos2_( 0 ) - pos1_( 0 ) ) +
                                std::atan2( pos2( 1 ) - pos1( 1 ), pos2( 0 ) - pos1( 0 ) ) );
        }

        return {out, angles};
    }

    std::pair<Eigen::MatrixX2d, std::vector<double>> cutTutteEmbedding( const topology::CutCombinatorialMap& map,
                                        const VertexPositionsFunc& vert_positions,
                                        const std::vector<std::pair<topology::Vertex, topology::Vertex>>& cut_extremities,
                                        const Laplace2dEdgeWeights& edge_weights_type )
    {
        if( map.dim() != 2 ) throw std::runtime_error( "Tutte embedding only supports 2d maps" );

        const auto edge_weights = [&]( const topology::Edge& e ) -> double {
            return edgeWeightLaplace2d( map, vert_positions, edge_weights_type, e );
        };

        const topology::CombinatorialMapBoundary bdry( map );
        const size_t n_bdry_verts = cellCount( bdry, 0 );
        const auto bdry_vertex_ids = indexingOrError( bdry, 0 );

        /*
           - For all the interior vertices of map, use the normal laplaceOperatorRowSparse.
           - For the cut vertices, use a modified laplaceOperatorRowSparse that
             acts as if both copies of the vertices and their neighborhood are at the same place.
             Also add another equation that constrains them to be the correct rotations of each other.
           - For the boundary (noncut) vertices, use constraints.
        */

        t.start( 0 );

        using SparseVectorXd = Eigen::SparseVector<double>;
        using SparseMatrixXd = Eigen::SparseMatrix<double>;
        std::map<size_t, Eigen::Index> unknown_verts;

        const size_t n_verts = cellCount( map, 0 );

        const auto vertex_ids = indexingOrError( map, 0 );
        const topology::CombinatorialMap& uncut_map = map.baseMap();
        const auto uncut_vertex_ids = indexingOrError( uncut_map, 0 );

        const auto other_side_of_cut = [&]( const topology::Vertex& v ) -> std::optional<topology::Vertex> {
            const topology::Vertex v_bdry = [&](){
                const auto maybe = maybeBoundaryDart( map, v );
                if( not maybe.has_value() )
                {
                    std::cerr << "Vertex " << v << " is not on the boundary!" << std::endl;
                    std::cerr << vert_positions( v ).transpose() << std::endl;
                }
                return maybe.value();
            }();

            return phi( uncut_map, {2,1}, v_bdry.dart() ).transform( [&]( const topology::Dart d ){
                return topology::Vertex( d );
            } );
        };

        const std::map<size_t, Eigen::Index> constrained_verts = [&]() -> std::map<size_t, Eigen::Index> {
            std::map<size_t, Eigen::Index> constrained;
            iterateCellsWhile( bdry, 0, [&]( const topology::Vertex& v ) {
                const topology::Vertex uncut_v = bdry.toUnderlyingCell( v );
                if( boundaryAdjacent( uncut_map, uncut_v ) )
                {
                    constrained.emplace( vertex_ids( uncut_v ), constrained.size() );
                }
                return true;
            } );
            return constrained;
        }();

        const size_t n_constrained_verts = constrained_verts.size();

        std::vector<Eigen::Triplet<double>> L_triplets;
        L_triplets.reserve( 4 * cellCount( map, 1 ) + 2 * n_verts + 5 * ( n_bdry_verts - n_constrained_verts ) );//TODO: update this

        const auto [BCs, cut_angles] = cutTutteBCsAndAngles( map, vert_positions, constrained_verts, cut_extremities );
        LOG( LOG_LAPLACE ) << "Cut angles: " << cut_angles << std::endl;

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
                                                     const Eigen::Index ref_vert_i_col,
                                                     const Eigen::Index ref_vert_j_col,
                                                     std::vector<Eigen::Triplet<double>>& triplets_out ) {
            LOG( LOG_LAPLACE ) << "Adding constrain rotation rows between " << i_col << " and " << j_col
                               << " with ref vert cols " << ref_vert_i_col << ", " << ref_vert_j_col << " at row "
                               << row << std::endl;
            triplets_out.emplace_back( row, 2 * i_col, -1 );
            triplets_out.emplace_back( row + 1, 2 * i_col + 1, -1 );

            triplets_out.emplace_back( row, 2 * j_col, rot.coeff( 0, 0 ) );
            triplets_out.emplace_back( row, 2 * j_col + 1, rot.coeff( 0, 1 ) );
            triplets_out.emplace_back( row + 1, 2 * j_col, rot.coeff( 1, 0 ) );
            triplets_out.emplace_back( row + 1, 2 * j_col + 1, rot.coeff( 1, 1 ) );

            triplets_out.emplace_back( row, 2 * ref_vert_i_col, 1 );
            triplets_out.emplace_back( row + 1, 2 * ref_vert_i_col + 1, 1 );

            triplets_out.emplace_back( row, 2 * ref_vert_j_col, - rot.coeff( 0, 0 ) );
            triplets_out.emplace_back( row, 2 * ref_vert_j_col + 1, -rot.coeff( 0, 1 ) );
            triplets_out.emplace_back( row + 1, 2 * ref_vert_j_col, -rot.coeff( 1, 0 ) );
            triplets_out.emplace_back( row + 1, 2 * ref_vert_j_col + 1, - rot.coeff( 1, 1 ) );
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
        for( size_t cut_ii = 0; cut_ii < cut_extremities.size(); cut_ii++ )
        {
            const auto& [start_v, end_v] = cut_extremities.at( cut_ii );
            // We take phi1s, add the rows for the vertices on both sides of the boundary, repeat until
            // we get to the next cone vertex, skip that, and continue to the last cone vertex.

            // Since we're treating these as boundary cmap darts, we don't actually process the cut starts, just the next vertex.
            topology::Dart curr_d = start_v.dart();
            SparseMatrixXd rot = sparse_rotation( cut_angles.at( cut_ii ) );

            const size_t comparison_vertex_column_i = vertex_ids( start_v );
            const size_t comparison_vertex_column_j = vertex_ids( other_side_of_cut( start_v ).value() );

            while( true )
            {
                const topology::Vertex curr_v( bdry.toUnderlyingCell( topology::Vertex( curr_d ) ) );
                if( boundaryAdjacent( uncut_map, curr_v ) ) break;

                const auto maybe_other_v = other_side_of_cut( curr_v );
                if( not maybe_other_v.has_value() ) throw std::runtime_error( "Expected other side of cut to exist" );
                const topology::Vertex other_v = maybe_other_v.value();
                const Eigen::Index i = unknown_verts.size();

                t.start( 2 );
                const SparseVectorXd row_1 = laplaceOperatorRowSparse( map, curr_v, edge_weights, n_verts );
                const SparseVectorXd row_2 = laplaceOperatorRowSparse( map, other_v, edge_weights, n_verts );
                t.stop( 2 );

                add_doubled_row( row_1, i, L_triplets );
                add_rotated_doubled_row( row_2, rot, i, L_triplets );
                add_constrain_rotation_rows( rot, 2 * i + 2, vertex_ids( curr_v ), vertex_ids( other_v ), comparison_vertex_column_i, comparison_vertex_column_j, L_triplets );

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

        return {result, cut_angles};
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

    std::map<size_t, Eigen::Vector2d>
        boundaryConstraints( const topology::CombinatorialMap& cut_cmap,
                             const VertexPositionsFunc& cut_cmap_positions,
                             const size_t n_cuts,
                             const std::function<bool( const topology::Vertex& )>& is_cut_extremity,
                             const std::function<bool( const topology::Vertex& )>& is_start_v )
    {
        if( cut_cmap.dim() != 2 ) throw std::invalid_argument( "boundaryConstraints only supports 2d cmaps" );
        const topology::CombinatorialMapBoundary bdry( cut_cmap );
        const auto bdry_positions = boundaryVertexPositions( bdry, cut_cmap_positions );
        const auto bdry_vert_ids = indexingOrError( bdry, 0 );
        const auto cut_vert_ids = indexingOrError( cut_cmap, 0 );

        const topology::Dart start_d = [&]() {
            std::optional<topology::Dart> d;
            iterateDartsWhile( bdry, [&]( const topology::Dart& a ) {
                if( is_start_v( bdry.toUnderlyingCell( topology::Vertex( a ) ) ) )
                {
                    d.emplace( a );
                    return false;
                }
                return true;
            } );
            if( not d.has_value() ) throw std::invalid_argument( "No vertex found for which is_start_v is true" );
            return d.value();
        }();

        std::map<size_t, Eigen::Vector2d> out;

        topology::Dart d = start_d;

        const size_t n_sides = n_cuts * 4;

        const std::vector<Eigen::Vector2d> ngon_verts = util::regularNGonVertices( n_sides );

        for( size_t side_ii = 0; side_ii < n_sides; side_ii++ )
        {
            std::map<topology::Vertex, double> side_positions;
            double cumulative_length = 0.0;
            do
            {
                side_positions.insert( { bdry.toUnderlyingCell( topology::Vertex( d ) ), cumulative_length } );
                cumulative_length += edgeLength( bdry, bdry_positions, d );
                d = phi( bdry, 1, d ).value();
            } while( not is_cut_extremity( bdry.toUnderlyingCell( topology::Vertex( d ) ) ) );

            const double factor = 1 / cumulative_length;
            for( auto& pr : side_positions )
            {
                const double s = pr.second * factor;
                out.emplace( cut_vert_ids( pr.first ), ( 1.0 - s ) * ngon_verts.at( side_ii ) + s * ngon_verts.at( side_ii + 1 ) );
            }
        }

        return out;
    }
} // namespace reparam