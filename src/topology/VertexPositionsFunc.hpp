#pragma once
#include<CustomEigen.hpp>

namespace topology
{
    class Vertex;
}

using VertexPositionsFunc = std::function<Vector3dMax( const topology::Vertex& )>;