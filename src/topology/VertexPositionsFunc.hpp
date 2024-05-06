#pragma once
#include<Eigen/Dense>

namespace topology
{
    class Vertex;
}

using VertexPositionsFunc = std::function<Eigen::Vector3d( const topology::Vertex& )>;