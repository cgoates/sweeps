#pragma once
#include <string>
#include <Eigen/Dense>

namespace cgogn
{
    struct CMap3;
}

namespace io
{
    void outputSimplicialFieldToVTK( const cgogn::CMap3& map,
                                     const Eigen::MatrixXd& data,
                                     const std::string& filename );
}