#pragma once
#include <Eigen/Dense>

using Vector3dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3>;
using Vector6dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 6>;
using MatrixX3dMax = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, 3>;