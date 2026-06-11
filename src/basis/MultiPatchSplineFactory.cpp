#include <MultiPatchSplineFactory.hpp>
#include <MultiPatchParametricAtlas.hpp>

namespace basis
{
    namespace
    {
        std::shared_ptr<const TPSplineSpace> buildTensorProductSplineSpace(
            const std::shared_ptr<const param::TPParametricAtlas>& atlas,
            const SmallVector<std::shared_ptr<const BSplineSpace1d>, 3>& bases_1d )
        {
            if( bases_1d.size() == 2 )
            {
                const auto bc = std::make_shared<const TPBasisComplex>(
                    atlas, bases_1d.at( 0 )->basisComplexPtr(), bases_1d.at( 1 )->basisComplexPtr() );
                return std::make_shared<const TPSplineSpace>( bc, bases_1d.at( 0 ), bases_1d.at( 1 ) );
            }

            if( bases_1d.size() == 3 )
            {
                const auto source_atlas =
                    std::dynamic_pointer_cast<const param::TPParametricAtlas>( atlas->sourcePtr() );
                if( not source_atlas )
                    throw std::runtime_error( "Cannot rebuild a 3D tensor-product spline space without a TP source atlas." );

                const auto source_bc = std::make_shared<const TPBasisComplex>(
                    source_atlas, bases_1d.at( 0 )->basisComplexPtr(), bases_1d.at( 1 )->basisComplexPtr() );
                const auto source_basis =
                    std::make_shared<const TPSplineSpace>( source_bc, bases_1d.at( 0 ), bases_1d.at( 1 ) );

                const auto bc =
                    std::make_shared<const TPBasisComplex>( atlas, source_bc, bases_1d.at( 2 )->basisComplexPtr() );
                return std::make_shared<const TPSplineSpace>( bc, source_basis, bases_1d.at( 2 ) );
            }

            throw std::runtime_error( "Only 2D and 3D tensor-product spline spaces are supported." );
        }
    }

    MultiPatchSplineSpace buildH1MultiPatchSplineSpace(
        const std::vector<std::shared_ptr<const TPSplineSpace>>& patches,
        const topology::MultiPatchCombinatorialMap::InternalConnectionsMap& connections )
    {
        return buildMultiPatchSplineSpace( patches, connections );
    }

    VectorConformingMultiPatchSplineSpace
        buildVectorConformingMultiPatchSplineSpace( const MultiPatchSplineSpace& primal,
                                                    const ConformingType conforming_type )
    {
        std::vector<std::shared_ptr<const VectorConformingTPSplineSpace>> vector_patches;
        vector_patches.reserve( primal.subSpaces().size() );

        for( const auto& patch : primal.subSpaces() )
        {
            const auto patch_bc =
                std::make_shared<const VectorConformingBasisComplex>( patch->basisComplexPtr(), conforming_type );
            vector_patches.push_back( std::make_shared<const VectorConformingTPSplineSpace>( patch_bc, *patch ) );
        }

        const auto bc =
            std::make_shared<const VectorConformingBasisComplex>( primal.basisComplexPtr(), conforming_type );
        return VectorConformingMultiPatchSplineSpace( bc, vector_patches );
    }

    VectorConformingMultiPatchSplineSpace buildHDivMultiPatchSplineSpace( const MultiPatchSplineSpace& primal )
    {
        return buildVectorConformingMultiPatchSplineSpace( primal, ConformingType::Divergence );
    }

    VectorConformingMultiPatchSplineSpace buildHCurlMultiPatchSplineSpace( const MultiPatchSplineSpace& primal )
    {
        return buildVectorConformingMultiPatchSplineSpace( primal, ConformingType::Curl );
    }

    MultiPatchSplineSpace buildL2MultiPatchSplineSpace( const VectorConformingMultiPatchSplineSpace& vector_space )
    {
        const auto& mp_atlas =
            dynamic_cast<const param::MultiPatchParametricAtlas&>( vector_space.basisComplex().parametricAtlas() );

        std::vector<std::shared_ptr<const TPSplineSpace>> l2_patches;
        l2_patches.reserve( vector_space.subSpaces().size() );

        for( const auto& vector_patch : vector_space.subSpaces() )
        {
            const auto& scalar_components = vector_patch->scalarTPBases();
            if( scalar_components.empty() )
                throw std::runtime_error( "Cannot build L2 multipatch space from a vector patch with no scalar components." );

            l2_patches.push_back( buildTensorProductSplineSpace(
                scalar_components.front()->basisComplex().parametricAtlasPtr(),
                vector_patch->reducedDegree1dBases() ) );
        }

        return buildDiscontinuousMultiPatchSplineSpace( l2_patches, mp_atlas.cmap().connections() );
    }
}
