#pragma once
#include <iostream>
#include <SimplicialComplex.hpp>
#include <VTKOutput.hpp>
#include<set>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>
#include <SweepInput.hpp>
#include<Logging.hpp>

#define LOG_LAPLACE 0

Eigen::Vector3d triangleNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3 )
{
    return ( v2 - v1 ).cross( v3 - v1 ).normalized();
}

Eigen::Vector3d triangleNormal( const cgogn::CMap3& map, const cgogn::CMap3::Face& f )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = f.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    const Eigen::Vector3d& pos3 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
    return triangleNormal( pos1, pos2, pos3 );
}

double edgeLength( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
{
    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    const cgogn::Dart& d = e.dart_;
    const Eigen::Vector3d& pos1 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( d ) );
    const Eigen::Vector3d& pos2 = cgogn::value<Eigen::Vector3d>( map, position, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
    return ( pos2 - pos1 ).norm();
}

double dihedralCotangent( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
{
    const Eigen::Vector3d n1 = triangleNormal( map, cgogn::CMap3::Face( e.dart_ ) );
    const Eigen::Vector3d n2 = triangleNormal( map, cgogn::CMap3::Face( cgogn::phi2( map, e.dart_ ) ) );
    //LOG( LOG_LAPLACE ) << "Dihedral angle between " << n1.transpose() << " and " << n2.transpose() << std::endl;

    const double cos_theta = abs( n1.dot( n2 ) );
    return cos_theta / std::sqrt( 1 - cos_theta * cos_theta );
}

double edgeWeight( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e )
{
    double weight = 0;
    const double factor = edgeLength( map, e ) / 12;
    //LOG( LOG_LAPLACE ) << "Edge " << e << " Factor: " << factor << std::endl;
    cgogn::foreach_incident_volume( map, e, [&]( cgogn::CMap3::Volume v ){
        weight += factor * dihedralCotangent( map, cgogn::CMap3::Edge( v.dart_ ) );
        // weight += factor * dihedralCotangent( map, cgogn::CMap3::Edge( cgogn::phi<1,2,-1>(map, v.dart_) ) );
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

Eigen::VectorXd solvelaplace( const cgogn::CMap3& map, const std::set<VertexId>& zero_bcs, const std::set<VertexId>& one_bcs )
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
            L_top.row( interior_verts.size() ) = laplaceOperatorRow( map, v );// <- This guy is getting too full...
            interior_verts.push_back( vid.id() );
        }
        return true;
    } );

    std::vector<Eigen::Index> all_verts( interior_verts.begin(), interior_verts.end() );
    for( const auto& a : boundary_verts ) all_verts.push_back( a );
    LOG( LOG_LAPLACE ) << all_verts << std::endl;
    LOG( LOG_LAPLACE ) << std::endl;

    LOG( LOG_LAPLACE ) << "L_top_prev:\n" << L_top << std::endl;

    // Permute to match reordered vertex indexing
    // TODO: Find a way to make this more explicit and readable?
    //L_top = L_top( Eigen::indexing::all, all_verts ).eval();

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

void mapFromInput( const SweepInput& sweep_input, cgogn::CMap3& map )
{
    cgogn::io::VolumeImportData import;
    import.reserve( sweep_input.points.size(), sweep_input.simplices.size() );

    for( const auto& tet : sweep_input.simplices )
    {
        import.volumes_types_.push_back( cgogn::io::VolumeType::Tetra );
        for( size_t i = 0; i < 4; i++ )
        {
            import.volumes_vertex_indices_.push_back( tet.vertex( i ).id() );
        }
    }

    import.vertex_position_ = sweep_input.points;

    import_volume_data( map, import );
}