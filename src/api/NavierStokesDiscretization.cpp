#include <NavierStokesDiscretization.hpp>
#include <KnotVector.hpp>

#include <IndexOperations.hpp>
#include <HierarchicalTPSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>

namespace api
{
    basis::DivConfTPSplineSpace buildHDIV( const basis::TPSplineSpace& H1 )
    {
        const auto HDIV_bc = std::make_shared<basis::DivConfBasisComplex>( H1.basisComplexPtr() );
        return basis::DivConfTPSplineSpace( HDIV_bc, H1 );
    }

    basis::TPSplineSpace buildL2( const std::shared_ptr<const param::TPParametricAtlas>& param,
                                  const basis::DivConfTPSplineSpace& HDIV )
    {
        const auto L2_bc =
            std::make_shared<basis::TPBasisComplex>( param,
                                                     HDIV.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
                                                     HDIV.reducedDegree1dBases().at( 1 )->basisComplexPtr() );
        return basis::TPSplineSpace( L2_bc, HDIV.reducedDegree1dBases().at( 0 ), HDIV.reducedDegree1dBases().at( 1 ) );
    }

    NavierStokesTPDiscretization::NavierStokesTPDiscretization( const basis::KnotVector& kv_s,
                                                                const basis::KnotVector& kv_t,
                                                                const size_t degree_s,
                                                                const size_t degree_t,
                                                                const Eigen::Matrix2Xd& cpts )
        : H1_ss( basis::buildBSpline( { kv_s, kv_t }, { degree_s, degree_t } ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap(), { topology::Dart( 0 ) } ),
          cpts( cpts ),
          H1( H1_ss, 2 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}

    std::vector<std::vector<topology::Cell>>
        leafElements( const std::vector<std::shared_ptr<const topology::TPCombinatorialMap>>& refinement_levels,
                      const std::vector<std::vector<std::pair<size_t, size_t>>>& elem_indices_to_refine )
    {
        const size_t dim = refinement_levels.front()->dim();

        std::vector<std::vector<topology::Cell>> out( refinement_levels.size(), std::vector<topology::Cell>() );

        // Start by adding all the first level elements to the output
        iterateCellsWhile( *refinement_levels.at( 0 ), refinement_levels.at( 0 )->dim(), [&]( const topology::Cell& c ) {
            out.front().push_back( c );
            return true;
        } );

        if( elem_indices_to_refine.size() == 0 ) return out;

        // Helper for iterating children
        const topology::HierarchicalTPCombinatorialMap temp_cmap( refinement_levels, out );

        // Convert the input indices to elements
        const auto elem_of_indices = [&refinement_levels,&dim]( const size_t level_ii, const std::pair<size_t, size_t>& iv ) {
            const SmallVector<topology::Dart, 3> unflat_darts( { topology::Dart( iv.first ), topology::Dart( iv.second ) } );
            return topology::Cell( flattenFull( *refinement_levels.at( level_ii ), topology::FullyUnflattenedDart( unflat_darts ) ), dim );
        };

        // Make changes according to the refinement for each level
        for( size_t level_ii = 0; level_ii < elem_indices_to_refine.size() - 1; level_ii++ )
        {
            auto& level_out = out.at( level_ii );
            const auto& level_refine = elem_indices_to_refine.at( level_ii );
            for( const std::pair<size_t, size_t>& elem_indices : level_refine )
            {
                const topology::Cell elem = elem_of_indices( level_ii, elem_indices );
                const auto it = std::remove( level_out.begin(), level_out.end(), elem );
                if( it == level_out.end() ) throw std::runtime_error( "Asked to refine an element that isn't active" );
                level_out.erase( it, level_out.end() );
                temp_cmap.iterateChildren( elem, level_ii, [&]( const topology::Cell& child_elem ) {
                    out.at( level_ii + 1 ).push_back( child_elem );
                    return true;
                } );
            }
        }

        return out;
    }

    basis::HierarchicalTPSplineSpace buildHierarchicalH1(
        const basis::KnotVector& kv_s,
        const basis::KnotVector& kv_t,
        const size_t degree_s,
        const size_t degree_t,
        const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine )
    {
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> refinement_levels;
        refinement_levels.reserve( elems_to_refine.size() + 1 );

        SmallVector<basis::KnotVector, 3> kvs( { kv_s, kv_t } );
        const SmallVector<size_t, 3> degrees( { degree_s, degree_t } );
        refinement_levels.push_back( std::make_shared<const basis::TPSplineSpace>(
            basis::buildBSpline( kvs, degrees ) ) );

        for( size_t i = 0; i < elems_to_refine.size(); i++ )
        {
            kvs.at( 0 ) = basis::dyadicRefine( kvs.at( 0 ) );
            kvs.at( 1 ) = basis::dyadicRefine( kvs.at( 1 ) );
            refinement_levels.push_back( std::make_shared<const basis::TPSplineSpace>(
                basis::buildBSpline( kvs, degrees ) ) );
        }

        std::vector<std::shared_ptr<const topology::TPCombinatorialMap>> refinement_level_cmaps;
        refinement_level_cmaps.reserve( refinement_levels.size() );
        std::transform( refinement_levels.begin(),
                        refinement_levels.end(),
                        std::back_inserter( refinement_level_cmaps ),
                        []( const auto& ss ) { return ss->basisComplex().parametricAtlas().cmapPtr(); } );

        return basis::buildHierarchicalSplineSpace( refinement_levels, leafElements( refinement_level_cmaps, elems_to_refine ) );
    }

    basis::DivConfHierarchicalTPSplineSpace buildHDIV( const basis::HierarchicalTPSplineSpace& H1 )
    {
        const auto HDIV_bc = std::make_shared<basis::DivConfBasisComplex>( H1.basisComplexPtr() );
        return basis::DivConfHierarchicalTPSplineSpace( HDIV_bc, H1 );
    }

    basis::HierarchicalTPSplineSpace buildL2( const std::shared_ptr<const param::HierarchicalTPParametricAtlas>& param,
                                  const basis::DivConfHierarchicalTPSplineSpace& HDIV )
    {
        std::vector<std::shared_ptr<const basis::TPBasisComplex>> refinement_level_bcs;
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> refinement_level_ss;
        refinement_level_bcs.reserve( param->cmap().numLevels() );
        refinement_level_ss.reserve( param->cmap().numLevels() );

        const auto& param_levels = param->refinementLevels();

        const auto& scalar_basis_levels_s = HDIV.scalarBases().at( 0 )->refinementLevels();
        const auto& scalar_basis_levels_t = HDIV.scalarBases().at( 1 )->refinementLevels();

        for( size_t i = 0; i < param->cmap().numLevels(); i++ )
        {
            const auto& source_ss = std::dynamic_pointer_cast<const basis::BSplineSpace1d>( scalar_basis_levels_t.at( i )->sourcePtr() );
            const auto& line_ss = scalar_basis_levels_s.at( i )->linePtr();
            if( source_ss.get() == nullptr ) throw std::runtime_error( "Bad cast to 1d bspline in creating L2 basis" );

            refinement_level_bcs.push_back( std::make_shared<const basis::TPBasisComplex>(
                param_levels.at( i ), source_ss->basisComplexPtr(), line_ss->basisComplexPtr() ) );

            refinement_level_ss.push_back( std::make_shared<const basis::TPSplineSpace>( refinement_level_bcs.back(), source_ss, line_ss ) );
        }

        const auto L2_bc = std::make_shared<basis::HierarchicalTPBasisComplex>( param, refinement_level_bcs );
        return basis::HierarchicalTPSplineSpace( L2_bc, refinement_level_ss );
    }

    NavierStokesHierarchicalDiscretization::NavierStokesHierarchicalDiscretization(
        const basis::KnotVector& kv_s,
        const basis::KnotVector& kv_t,
        const size_t degree_s,
        const size_t degree_t,
        const Eigen::Matrix2Xd& unrefined_cpts,
        const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine )
        : H1_ss( buildHierarchicalH1( kv_s, kv_t, degree_s, degree_t, elems_to_refine ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap() ),
          cpts( unrefined_cpts * prolongationOperator( H1_ss ).transpose() ),
          H1( H1_ss, 2 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}
}