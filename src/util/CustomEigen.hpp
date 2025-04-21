#pragma once
#include <Eigen/Core>

using Vector1d = Eigen::Matrix<double, 1, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Vector3dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 3>;
using Vector6dMax = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 6>;
using MatrixX3dMax = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, Eigen::Dynamic, 3>;