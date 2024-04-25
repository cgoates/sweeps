#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Simplex.hpp>

struct SimplicialComplex
{
    std::vector<Simplex> simplices;
    std::vector<Eigen::Vector3d> points;
};
