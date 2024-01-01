#include<Laplace.hpp>
#include<Simplex.hpp>
#include<Logging.hpp>
#include <cgogn/core/types/cell_marker.h>
#include<SimplexUtilities.hpp>

#define LOG_LAPLACE 0

double edgeWeight( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
{
    double weight = 0;
    const double factor = SimplexUtilities::edgeLength( map, e ) / 12;
    LOG( LOG_LAPLACE ) << "Edge " << e << " Factor: " << factor << std::endl;
    cgogn::foreach_incident_volume( map, e, [&]( cgogn::CMap3::Volume v ){
        weight += factor * SimplexUtilities::dihedralCotangent( map, cgogn::CMap3::Edge( v.dart_ ) );
        // weight += factor * SimplexUtilities::dihedralCotangent( map, cgogn::CMap3::Edge( cgogn::phi<1,2,-1>(map, v.dart_) ) );
        return true;
    } );

    return weight;
}

// FIXME: This isn't right; I think I'm using the wrong dihedral angle?
Eigen::VectorXd laplaceOperatorRow( const cgogn::CMap3& map, const cgogn::CMap3::Vertex& v1 )
{
    const int n_verts = cgogn::nb_cells<cgogn::CMap3::Vertex>( map );
    LOG( LOG_LAPLACE ) << n_verts << std::endl;
    Eigen::VectorXd out = Eigen::VectorXd::Zero( n_verts );
    const VertexId vid1 = cgogn::index_of( map, v1 );
    cgogn::foreach_incident_edge( map, v1, [&]( cgogn::CMap3::Edge e ){
        const double edge_weight = edgeWeight( map, e );
        const VertexId vid2 = cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, e.dart_ ) ) );
        LOG( LOG_LAPLACE ) << vid1.id() << "->" << vid2.id() << " : " << edge_weight << std::endl;

        out( vid1.id() ) -= edge_weight;
        out( vid2.id() ) += edge_weight;
        return true;
    } );

    return out;
}

Eigen::MatrixXd laplaceOperator( const cgogn::CMap3& map )
{
    const int n_verts = cgogn::nb_cells<cgogn::CMap3::Vertex>( map );
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero( n_verts, n_verts );
    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Vertex v1 ) {
        const VertexId vid1 = cgogn::index_of( map, v1 );
        out.row( vid1.id() ) = laplaceOperatorRow( map, v1 );
        return true;
    } );

    return out;
}

Eigen::VectorXd solveLaplace( const cgogn::CMap3& map, const std::set<VertexId>& zero_bcs, const std::set<VertexId>& one_bcs )
{
    std::vector<Eigen::Index> interior_verts;
    std::vector<Eigen::Index> boundary_verts;

    const int n_verts = cgogn::nb_cells<cgogn::CMap3::Vertex>( map );

    Eigen::VectorXd BCs = Eigen::VectorXd::Zero( zero_bcs.size() + one_bcs.size() );
    Eigen::MatrixXd L_top( n_verts - BCs.size(), n_verts );

    LOG( LOG_LAPLACE ) << "About to fill rows. " << L_top.rows() << ", " << L_top.cols() << "\n";
    LOG( LOG_LAPLACE ) << "Also BCs. " << BCs.rows() << std::endl;
    LOG( LOG_LAPLACE ) << "zeros: " << zero_bcs << std::endl;
    LOG( LOG_LAPLACE ) << "ones: " << one_bcs << std::endl;
    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Vertex v ) {
        const VertexId vid = cgogn::index_of( map, v );
        if( zero_bcs.contains( vid ) )
        {
            LOG( LOG_LAPLACE ) << "zero bc\n";
            boundary_verts.push_back( vid.id() );
        }
        else if( one_bcs.contains( vid ) )
        {
            LOG( LOG_LAPLACE ) << "one bc\n";
            BCs( boundary_verts.size() ) = 1.0;
            boundary_verts.push_back( vid.id() );
        }
        else
        {
            LOG( LOG_LAPLACE ) << "interior: " << interior_verts.size() << "\n";
            LOG( LOG_LAPLACE ) << interior_verts << std::endl;
            L_top.row( interior_verts.size() ) = laplaceOperatorRow( map, v );
            interior_verts.push_back( vid.id() );
        }
        return true;
    } );

    LOG( LOG_LAPLACE ) << "L_top: " << std::endl << L_top << std::endl << std::endl;
    LOG( LOG_LAPLACE ) << "BCs: " << std::endl << BCs << std::endl << std::endl;

    const Eigen::MatrixXd L_II = L_top( Eigen::indexing::all, interior_verts );
    const Eigen::VectorXd rhs = -L_top( Eigen::indexing::all, boundary_verts ) * BCs;

    LOG( LOG_LAPLACE ) << "L_II:\n" << L_II << std::endl << std::endl;
    LOG( LOG_LAPLACE ) << "rhs:\n" << rhs << std::endl << std::endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver( L_II );
    const Eigen::VectorXd ans = solver.solve( rhs );

    Eigen::VectorXd result( L_top.cols() );
    result( interior_verts ) = ans;
    result( boundary_verts ) = BCs;

    return result;
}
