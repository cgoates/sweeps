#pragma once
#include <string>
#include <Eigen/Dense>

namespace topology
{
    class MultiPatchCombinatorialMap;
    class Cell;
}
namespace basis
{
    class MultiPatchSplineSpace;
}

namespace io
{
    void outputMultiPatchSplinesToMFEM( const basis::MultiPatchSplineSpace& mpss,
                                        const Eigen::Ref<const Eigen::MatrixXd> geom,
                                        const std::function<size_t( const size_t )>& patch_to_block_id,
                                        const std::function<size_t( const topology::MultiPatchCombinatorialMap&, const topology::Cell& )>& bdry_cell_to_block_id,
                                        const std::string& filename );

    void outputMultiPatchSplinesToMFEM( const basis::MultiPatchSplineSpace& mpss,
                                        const Eigen::Ref<const Eigen::MatrixXd> geom,
                                        const std::string& filename );
}