#pragma once
#include<Eigen/Dense>
#include<set>
#include <cgogn/core/types/maps/cmap/cmap3.h>

class VertexId;

double edgeWeight( const cgogn::CMap3& map, const cgogn::CMap3::Edge& e );

Eigen::VectorXd laplaceOperatorRow( const cgogn::CMap3& map, const cgogn::CMap3::Vertex& v1 );

Eigen::MatrixXd laplaceOperator( const cgogn::CMap3& map );

Eigen::VectorXd solveLaplace( const cgogn::CMap3& map, const std::set<VertexId>& zero_bcs, const std::set<VertexId>& one_bcs );

Eigen::VectorXd solveLaplaceSparse( const cgogn::CMap3& map, const std::set<VertexId>& zero_bcs, const std::set<VertexId>& one_bcs );