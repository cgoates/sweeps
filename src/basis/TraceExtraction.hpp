#pragma once
#include <ParentBasis.hpp>
#include <SplineSpace.hpp>
#include <TPCombinatorialMap.hpp>
#include <Eigen/Core>
#include <optional>

namespace basis
{
    /// Identifies one side of a tensor-product element.
    ///
    /// The side-id convention matches MultiPatchCombinatorialMap: side 2*axis is the
    /// upper/coordinate-1 side and side 2*axis+1 is the lower/coordinate-0 side.
    struct ElementSide
    {
        size_t axis = 0;
        bool lower = true;

        ElementSide() = default;
        ElementSide( const size_t axis_in, const bool lower_in ) : axis( axis_in ), lower( lower_in ) {}

        size_t sideId() const { return 2 * axis + ( lower ? 1 : 0 ); }
        ElementSide opposite() const { return ElementSide( axis, not lower ); }
        bool operator==( const ElementSide& ) const = default;
    };

    ElementSide elementSideFromId( const size_t side_id );

    struct TraceSideData
    {
        topology::Cell element;
        ElementSide side;
        std::vector<FunctionId> connectivity;
        Eigen::MatrixXd extraction;
        std::vector<size_t> element_row_indices;
    };

    struct TraceInterfaceElement
    {
        TraceSideData first;
        std::optional<TraceSideData> second;

        bool isBoundary() const { return not second.has_value(); }
    };

    SmallVector<size_t, 3> traceDegrees( const SmallVector<size_t, 3>& degrees, const ElementSide& side );

    std::vector<Eigen::Index> traceColumnIndices( const SmallVector<size_t, 3>& degrees,
                                                  const ElementSide& side );

    std::vector<Eigen::Index> traceColumnIndices( const ParentBasis& pb, const ElementSide& side );

    Eigen::MatrixXd traceExtractionMatrix( const Eigen::MatrixXd& element_extraction,
                                           const ParentBasis& pb,
                                           const ElementSide& side );

    TraceSideData traceSideData( const SplineSpace& ss,
                                 const topology::Cell& element,
                                 const ElementSide& side,
                                 const double row_tol = 1e-12 );

    std::optional<std::pair<topology::Cell, ElementSide>>
        adjacentElementSide( const topology::TPCombinatorialMap& cmap,
                             const topology::Cell& element,
                             const ElementSide& side );

    TraceInterfaceElement traceInterfaceElement( const SplineSpace& ss,
                                                 const topology::Cell& element,
                                                 const ElementSide& side,
                                                 const double row_tol = 1e-12 );

    TraceInterfaceElement traceInterfaceElement( const SplineSpace& first_ss,
                                                 const topology::Cell& first_element,
                                                 const ElementSide& first_side,
                                                 const SplineSpace& second_ss,
                                                 const topology::Cell& second_element,
                                                 const ElementSide& second_side,
                                                 const double row_tol = 1e-12 );
}
