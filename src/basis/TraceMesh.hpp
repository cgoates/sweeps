#pragma once
#include <TraceExtraction.hpp>
#include <MultiPatchCombinatorialMap.hpp>
#include <Eigen/Core>
#include <optional>
#include <vector>

namespace basis
{
    enum class TraceMeshEntityType
    {
        Boundary,
        PatchInterface,
        Interior
    };

    enum class InteriorTraceMode
    {
        All,
        Broken
    };

    struct TraceMeshSide
    {
        topology::Cell element;
        ElementSide side;
        size_t patch_id = 0;
        size_t level = 0;
    };

    struct TraceMeshInterface
    {
        TraceMeshEntityType type;
        TraceMeshSide first;
        std::optional<TraceMeshSide> second;
        std::optional<topology::MultiPatchCombinatorialMap::TPPermutation> permutation;

        bool isBoundary() const { return not second.has_value(); }
    };

    std::vector<size_t> permuteTraceSideIndex(
        const std::vector<size_t>& index,
        const std::vector<size_t>& lengths,
        topology::MultiPatchCombinatorialMap::TPPermutation permutation );

    Eigen::VectorXd permuteTraceSidePoint(
        const Eigen::VectorXd& point,
        topology::MultiPatchCombinatorialMap::TPPermutation permutation );

    std::vector<TraceMeshInterface> boundaryTraceMeshInterfaces( const SplineSpace& ss );

    std::vector<TraceMeshInterface> patchTraceMeshInterfaces( const SplineSpace& ss );

    /// Enumerates same-patch, same-level leaf interfaces.  In hierarchical spaces this
    /// intentionally skips coarse-fine subdivision; that belongs with the physical trace
    /// integration layer, where quadrature-side subdivision can be handled explicitly.
    std::vector<TraceMeshInterface> interiorTraceMeshInterfaces(
        const SplineSpace& ss,
        InteriorTraceMode mode = InteriorTraceMode::All,
        double row_tol = 1e-12 );

    std::vector<TraceMeshInterface> traceMeshInterfaces(
        const SplineSpace& ss,
        bool include_boundary,
        bool include_patch_interfaces,
        bool include_interior,
        InteriorTraceMode interior_mode = InteriorTraceMode::All,
        double row_tol = 1e-12 );

    TraceInterfaceElement traceInterfaceElement(
        const SplineSpace& ss,
        const TraceMeshInterface& iface,
        double row_tol = 1e-12 );
}
