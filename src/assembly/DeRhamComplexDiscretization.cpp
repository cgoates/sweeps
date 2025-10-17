#include <DeRhamComplexDiscretization.hpp>
#include <KnotVector.hpp>
#include <IndexOperations.hpp>
#include <HierarchicalTPSplineSpace.hpp>
#include <CombinatorialMapMethods.hpp>

namespace assembly
{
    basis::VectorConformingTPSplineSpace buildHDIV( const basis::TPSplineSpace& H1 )
    {
        const auto HDIV_bc = std::make_shared<basis::VectorConformingBasisComplex>( H1.basisComplexPtr() );
        return basis::VectorConformingTPSplineSpace( HDIV_bc, H1 );
    }

    basis::VectorConformingTPSplineSpace buildHCURL( const basis::TPSplineSpace& H1 )
    {
        const auto HCURL_bc = std::make_shared<basis::VectorConformingBasisComplex>( H1.basisComplexPtr(), basis::ConformingType::Curl );
        return basis::VectorConformingTPSplineSpace( HCURL_bc, H1, basis::ConformingType::Curl );
    }

    basis::TPSplineSpace buildL2( const std::shared_ptr<const param::TPParametricAtlas>& param,
                                  const basis::VectorConformingTPSplineSpace& HDIV )
    {
        const auto L2_bc =
            std::make_shared<basis::TPBasisComplex>( param,
                                                     HDIV.reducedDegree1dBases().at( 0 )->basisComplexPtr(),
                                                     HDIV.reducedDegree1dBases().at( 1 )->basisComplexPtr() );
        return basis::TPSplineSpace( L2_bc, HDIV.reducedDegree1dBases().at( 0 ), HDIV.reducedDegree1dBases().at( 1 ) );
    }

    DeRhamComplexTPDiscretization::DeRhamComplexTPDiscretization( const basis::KnotVector& kv_s,
                                                                const basis::KnotVector& kv_t,
                                                                const size_t degree_s,
                                                                const size_t degree_t,
                                                                const Eigen::Matrix2Xd& cpts )
        : H1_ss( basis::buildBSpline( { kv_s, kv_t }, { degree_s, degree_t } ) ),
          HCURL_ss( buildHCURL( H1_ss ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap(), { topology::Dart( 0 ) } ),
          cpts( cpts ),
          H1( H1_ss, 2 ),
          HCURL( HCURL_ss, 1 ),
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
        for( size_t level_ii = 0; level_ii < elem_indices_to_refine.size(); level_ii++ )
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

    std::vector<std::vector<topology::Cell>>
        leafElements( const std::vector<std::shared_ptr<const topology::TPCombinatorialMap>>& refinement_levels,
                      const std::vector<std::vector<std::array<size_t, 3>>>& elem_indices_to_refine )
    {
        const size_t dim = refinement_levels.front()->dim();

        std::vector<std::vector<topology::Cell>> out( refinement_levels.size(), std::vector<topology::Cell>() );

        // Start by adding all the first level elements to the output
        iterateCellsWhile( *refinement_levels.at( 0 ), dim, [&]( const topology::Cell& c ) {
            out.front().push_back( c );
            return true;
        } );

        if( elem_indices_to_refine.size() == 0 ) return out;

        // Helper for iterating children
        const topology::HierarchicalTPCombinatorialMap temp_cmap( refinement_levels, out );

        // Convert the input indices to elements
        const auto elem_of_indices = [&refinement_levels,&dim]( const size_t level_ii, const std::array<size_t, 3>& iv ) {
            const SmallVector<topology::Dart, 3> unflat_darts( { topology::Dart( iv.at( 0 ) ), topology::Dart( iv.at( 1 ) ), topology::Dart( iv.at( 2 ) ) } );
            return topology::Cell( flattenFull( *refinement_levels.at( level_ii ), topology::FullyUnflattenedDart( unflat_darts ) ), dim );
        };

        // Make changes according to the refinement for each level
        for( size_t level_ii = 0; level_ii < elem_indices_to_refine.size(); level_ii++ )
        {
            auto& level_out = out.at( level_ii );
            const auto& level_refine = elem_indices_to_refine.at( level_ii );
            for( const std::array<size_t, 3>& elem_indices : level_refine )
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

    // std::vector<std::vector<std::array<size_t, 3>>>
    //     enforceAdmissibility( const std::vector<std::shared_ptr<const topology::TPCombinatorialMap>>& refinement_levels,
    //                           std::vector<std::vector<std::array<size_t, 3>>> elem_indices_to_refine,
    //                           const size_t admissibility,
    //                           const size_t degree )
    // {
    //     const size_t dim = refinement_levels.front()->dim();

    //     if( elem_indices_to_refine.size() == 0 ) return elem_indices_to_refine;

    //     const auto get_tp_lengths = []( const topology::TPCombinatorialMap& cmap ) -> util::IndexVec {
    //         const auto cmaps = tensorProductComponentCMaps( cmap );
    //         util::IndexVec out;
    //         std::transform( cmaps.begin(), cmaps.end(), std::back_inserter( out ), []( const auto& c ) { return topology::cellCount( *c, 1 ); } );
    //         return out;
    //     };

    //     std::vector<util::IndexVec> lengths;
    //     std::transform( refinement_levels.begin(), refinement_levels.end(), std::back_inserter( lengths ), get_tp_lengths );

    //     // Enforce admissibility
    //     for( int level_ii = elem_indices_to_refine.size() - 1; level_ii >= 0; level_ii-- )
    //     {
    //         for( const std::array<size_t, 3>& elem_indices : elem_indices_to_refine.at( level_ii ) )
    //         {

    //             for( size_t dlevel = 1; dlevel <= admissibility; dlevel++ )
    //             {
    //                 if( level_ii - dlevel < 0 ) break;
    //                 std::array<size_t, 3> coarse_parent;
    //                 std::transform( elem_indices.begin(),
    //                                 elem_indices.end(),
    //                                 coarse_parent.begin(),
    //                                 [ratio = lengths.at( level_ii ).at( 0 ) / lengths.at( level_ii - dlevel ).at( 0 )](
    //                                     size_t x ) { return x / ratio; } );
    //                 const size_t radius = std::ceil( static_cast<double>( degree + 1 ) / std::pow( 2.0, dlevel ) );
    //                 auto& coarser_level_refine = elem_indices_to_refine.at( level_ii - dlevel );

    //                 util::iterateTensorProduct( util::IndexVec( dim, 2*radius + 1 ), [&]( const util::IndexVec& it ) {
    //                     // If out of bounds, skip it
    //                     bool out_of_bounds = false;
    //                     std::array<size_t, 3> coarse_neighbor;
    //                     for( size_t i = 0; i < dim; i++ )
    //                     {
    //                         coarse_neighbor.at( i ) = coarse_parent.at( i ) + it.at( i ) - radius;
    //                         if( coarse_neighbor.at( i ) > lengths.at( level_ii - dlevel ).at( i ) )
    //                             out_of_bounds = true;
    //                     }
    //                     if( not out_of_bounds )
    //                         coarser_level_refine.push_back( coarse_neighbor );
    //                 } );
    //             }
    //         }
    //     }

    //     // Remove duplicates
    //     for( auto& level_refine : elem_indices_to_refine )
    //     {
    //         std::sort( level_refine.begin(), level_refine.end() );
    //         level_refine.erase( std::unique( level_refine.begin(), level_refine.end() ), level_refine.end() );
    //     }

    //     return elem_indices_to_refine;
    // }

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

    basis::HierarchicalTPSplineSpace buildHierarchicalH1(
        const basis::KnotVector& kv_s,
        const basis::KnotVector& kv_t,
        const basis::KnotVector& kv_u,
        const size_t degree_s,
        const size_t degree_t,
        const size_t degree_u,
        const std::vector<std::vector<std::array<size_t, 3>>>& elems_to_refine )
    {
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> refinement_levels;
        refinement_levels.reserve( elems_to_refine.size() + 1 );

        SmallVector<basis::KnotVector, 3> kvs( { kv_s, kv_t, kv_u } );
        const SmallVector<size_t, 3> degrees( { degree_s, degree_t, degree_u } );
        refinement_levels.push_back( std::make_shared<const basis::TPSplineSpace>(
            basis::buildBSpline( kvs, degrees ) ) );

        for( size_t i = 0; i < elems_to_refine.size(); i++ )
        {
            kvs.at( 0 ) = basis::dyadicRefine( kvs.at( 0 ) );
            kvs.at( 1 ) = basis::dyadicRefine( kvs.at( 1 ) );
            kvs.at( 2 ) = basis::dyadicRefine( kvs.at( 2 ) );
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

    basis::VectorConformingHierarchicalTPSplineSpace buildHDIV( const basis::HierarchicalTPSplineSpace& H1 )
    {
        const auto HDIV_bc = std::make_shared<basis::VectorConformingBasisComplex>( H1.basisComplexPtr() );
        return basis::VectorConformingHierarchicalTPSplineSpace( HDIV_bc, H1 );
    }

    basis::VectorConformingHierarchicalTPSplineSpace buildHCURL( const basis::HierarchicalTPSplineSpace& H1 )
    {
        const auto HCURL_bc = std::make_shared<basis::VectorConformingBasisComplex>( H1.basisComplexPtr(), basis::ConformingType::Curl );
        return basis::VectorConformingHierarchicalTPSplineSpace( HCURL_bc, H1 );
    }

    basis::HierarchicalTPSplineSpace buildL2( const std::shared_ptr<const param::HierarchicalTPParametricAtlas>& param,
                                              const basis::VectorConformingHierarchicalTPSplineSpace& HDIV )
    {
        std::vector<std::shared_ptr<const basis::TPBasisComplex>> refinement_level_bcs;
        std::vector<std::shared_ptr<const basis::TPSplineSpace>> refinement_level_ss;
        refinement_level_bcs.reserve( param->cmap().numLevels() );
        refinement_level_ss.reserve( param->cmap().numLevels() );

        const auto& param_levels = param->refinementLevels();

        const auto& scalar_basis_levels_s = HDIV.scalarBases().at( 0 )->refinementLevels();
        const auto& scalar_basis_levels_t = HDIV.scalarBases().at( 1 )->refinementLevels();

        const size_t dim = param->cmap().dim();

        if( dim == 2 )
        {
            for( size_t i = 0; i < param->cmap().numLevels(); i++ )
            {
                const auto& source_ss = std::dynamic_pointer_cast<const basis::BSplineSpace1d>(
                    scalar_basis_levels_t.at( i )->sourcePtr() );
                const auto& line_ss = scalar_basis_levels_s.at( i )->linePtr();
                if( source_ss.get() == nullptr )
                    throw std::runtime_error( "Bad cast to 1d bspline in creating L2 basis" );

                refinement_level_bcs.push_back( std::make_shared<const basis::TPBasisComplex>(
                    param_levels.at( i ), source_ss->basisComplexPtr(), line_ss->basisComplexPtr() ) );

                refinement_level_ss.push_back(
                    std::make_shared<const basis::TPSplineSpace>( refinement_level_bcs.back(), source_ss, line_ss ) );
            }
        }
        else if( dim == 3 )
        {
            for( size_t i = 0; i < param->cmap().numLevels(); i++ )
            {
                const auto& s_1d = std::dynamic_pointer_cast<const basis::BSplineSpace1d>(
                    std::dynamic_pointer_cast<const basis::TPSplineSpace>( scalar_basis_levels_t.at( i )->sourcePtr() )
                        ->sourcePtr() );
                const auto& t_1d =
                    std::dynamic_pointer_cast<const basis::TPSplineSpace>( scalar_basis_levels_s.at( i )->sourcePtr() )
                        ->linePtr();
                const auto& u_1d = scalar_basis_levels_s.at( i )->linePtr();
                if( s_1d.get() == nullptr or t_1d.get() == nullptr )
                    throw std::runtime_error( "Bad cast to 1d bspline in creating L2 basis" );

                const auto source_bc = std::make_shared<const basis::TPBasisComplex>(
                    std::dynamic_pointer_cast<const param::TPParametricAtlas>( param_levels.at( i )->sourcePtr() ),
                    s_1d->basisComplexPtr(),
                    t_1d->basisComplexPtr() );

                refinement_level_bcs.push_back( std::make_shared<const basis::TPBasisComplex>(
                    param_levels.at( i ), source_bc, u_1d->basisComplexPtr() ) );

                const auto source_ss = std::make_shared<const basis::TPSplineSpace>( source_bc, s_1d, t_1d );

                refinement_level_ss.push_back(
                    std::make_shared<const basis::TPSplineSpace>( refinement_level_bcs.back(), source_ss, u_1d ) );
            }
        }

        const auto L2_bc = std::make_shared<basis::HierarchicalTPBasisComplex>( param, refinement_level_bcs );
        return basis::HierarchicalTPSplineSpace( L2_bc, refinement_level_ss );
    }

    DeRhamComplexHierarchicalDiscretization::DeRhamComplexHierarchicalDiscretization(
        const basis::KnotVector& kv_s,
        const basis::KnotVector& kv_t,
        const size_t degree_s,
        const size_t degree_t,
        const Eigen::Matrix2Xd& unrefined_cpts,
        const std::vector<std::vector<std::pair<size_t, size_t>>>& elems_to_refine )
        : H1_ss( buildHierarchicalH1( kv_s, kv_t, degree_s, degree_t, elems_to_refine ) ),
          HCURL_ss( buildHCURL( H1_ss ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap() ),
          cpts( unrefined_cpts * prolongationOperator( H1_ss ).transpose() ),
          H1( H1_ss, 2 ),
          HCURL( HCURL_ss, 1 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}

    DeRhamComplexHierarchicalDiscretization::DeRhamComplexHierarchicalDiscretization(
        const basis::KnotVector& kv_s,
        const basis::KnotVector& kv_t,
        const basis::KnotVector& kv_u,
        const size_t degree_s,
        const size_t degree_t,
        const size_t degree_u,
        const Eigen::Matrix3Xd& unrefined_cpts,
        const std::vector<std::vector<std::array<size_t, 3>>>& elems_to_refine )
        : H1_ss( buildHierarchicalH1( kv_s, kv_t, kv_u, degree_s, degree_t, degree_u, elems_to_refine ) ),
          HCURL_ss( buildHCURL( H1_ss ) ),
          HDIV_ss( buildHDIV( H1_ss ) ),
          L2_ss( buildL2( H1_ss.basisComplex().parametricAtlasPtr(), HDIV_ss ) ),
          cmap_bdry( H1_ss.basisComplex().parametricAtlas().cmap() ),
          cpts( unrefined_cpts * prolongationOperator( H1_ss ).transpose() ),
          H1( H1_ss, 2 ),
          HCURL( HCURL_ss, 1 ),
          HDIV( HDIV_ss, 1 ),
          L2( L2_ss, 1 )
    {}

    void localizeElement( DeRhamComplexDiscretization& nsd, const topology::Cell& elem )
    {
        nsd.getH1().localizeElement( elem );
        nsd.getHCURL().localizeElement( elem );
        nsd.getHDIV().localizeElement( elem );
        nsd.getL2().localizeElement( elem );
    }

    void localizePoint( DeRhamComplexDiscretization& nsd, const param::ParentPoint& ppt )
    {
        nsd.getH1().localizePoint( ppt );
        nsd.getHCURL().localizePoint( ppt );
        nsd.getHDIV().localizePoint( ppt );
        nsd.getL2().localizePoint( ppt );
    }
}