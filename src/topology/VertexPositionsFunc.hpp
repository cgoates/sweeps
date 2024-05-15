#pragma once
#include<Eigen/Dense>

namespace topology
{
    class Vertex;
}

using Vector3dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3>;
using VertexPositionsFunc = std::function<Vector3dMax( const topology::Vertex& )>;