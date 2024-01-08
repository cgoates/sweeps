#pragma once
#include <Eigen/Dense>
#include <set>

class VertexId;
namespace cgogn
{
    struct CMap3;
}

Eigen::VectorXd solveLaplaceSparse( const cgogn::CMap3& map,
                                    const std::set<VertexId>& zero_bcs,
                                    const std::set<VertexId>& one_bcs );