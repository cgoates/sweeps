#pragma once
#include <MultiPatchSplineSpace.hpp>
#include <VectorConformingMultiPatchSplineSpace.hpp>

namespace basis
{
    MultiPatchSplineSpace buildH1MultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const topology::MultiPatchCombinatorialMap::InternalConnectionsMap& connections );

    VectorConformingMultiPatchSplineSpace
        buildVectorConformingMultiPatchSplineSpace( const MultiPatchSplineSpace& primal,
                                                    const ConformingType conforming_type );

    VectorConformingMultiPatchSplineSpace buildHDivMultiPatchSplineSpace( const MultiPatchSplineSpace& primal );

    VectorConformingMultiPatchSplineSpace buildHCurlMultiPatchSplineSpace( const MultiPatchSplineSpace& primal );

    MultiPatchSplineSpace buildL2MultiPatchSplineSpace( const VectorConformingMultiPatchSplineSpace& vector_space );
}
