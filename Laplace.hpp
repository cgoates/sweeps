#pragma once
#include <Eigen/Dense>
#include <set>

class VertexId;
class Normal;

namespace cgogn
{
    struct CMap3;
}

Eigen::VectorXd solveLaplaceSparse( const cgogn::CMap3& map,
                                    const std::vector<bool>& zero_bcs,
                                    const std::vector<bool>& one_bcs,
                                    const std::vector<Normal>& normals );