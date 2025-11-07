#include <ParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>
#include <CommonUtils.hpp>

namespace param
{
    BaryCoordIsZeroVec parentDomainBoundary( const ParametricAtlas& atlas, const topology::Cell& cell )
    {
        if( cell.dim() == 0 ) return atlas.parentPoint( cell ).mBaryCoordIsZero;
        else
        {
            BaryCoordIsZeroVec out( numTotalCoordinates( atlas.parentDomain( topology::Cell( cell.dart(), atlas.cmap().dim() ) ) ), true );
            iterateAdjacentCellsOfRestrictedCell( atlas.cmap(), cell, atlas.cmap().dim(), 0, [&]( const topology::Vertex& v ) {
                const ParentPoint ppt = atlas.parentPoint( v );
                out = join( out, ppt.mBaryCoordIsZero );
                return true;
            } );
            return out;
        }
    }

    size_t parametricLengthIndexAlongEdge( const ParametricAtlas& atlas, const topology::Edge& e )
    {
        const ParentDomain pd = atlas.parentDomain( topology::Cell( e.dart(), atlas.cmap().dim() ) );
        const BaryCoordIsZeroVec edge_bdry = parentDomainBoundary( atlas, e );

        size_t index = 0;
        SmallVector<size_t, 2> group_local_changing_indices;
        const size_t param_group_containing_edge = [&]() -> size_t {
            if( pd.coordinateGroups().size() == 1 ) return 0;

            size_t first_expanded_coord = 0;
            for( size_t group_ii = 0; group_ii < pd.coordinateGroups().size(); group_ii++ )
            {
                const CoordinateSystem& cs = pd.coordinateGroups().at( group_ii );

                SmallVector<size_t, 2> changing_indices;

                size_t n_nonzeros = 0;
                for( size_t coord = first_expanded_coord, e = first_expanded_coord + cs.dim(); coord <= e; coord++ )
                {
                    if( not edge_bdry.at( coord ) )
                    {
                        n_nonzeros++;
                        changing_indices.push_back( coord - first_expanded_coord );
                    }
                }

                if( n_nonzeros == 2 )
                {
                    group_local_changing_indices = changing_indices;
                    return group_ii;
                }
                else index += numParametricLengths( cs );

                first_expanded_coord += numTotalCoordinates( cs );
            }

            throw std::runtime_error( "No group found that contains the edge" );
        }();

        const size_t dim = pd.coordinateGroups().at( param_group_containing_edge ).dim();

        if( dim == 1 ) return index;

        const size_t k0 = group_local_changing_indices.at( 0 );
        const size_t k1 = group_local_changing_indices.at( 1 );

        //Given the two vertex indices in the group, calculates the lex ordering position.
        index += k0 * ( 2 * dim - k0 + 1 ) / 2 + k1 - k0 - 1;

        return index;
    }

    topology::Vertex originVertex( const ParametricAtlas& atlas, const topology::Cell& c )
    {
        if( c.dim() != atlas.cmap().dim() ) throw std::runtime_error( "Cell dimension does not match atlas dimension" );
        const ParentDomain pd = atlas.parentDomain( c );

        std::optional<topology::Vertex> out;
        iterateAdjacentCellsOfRestrictedCell( atlas.cmap(), c, atlas.cmap().dim(), 0, [&]( const topology::Vertex& v ) {
            const param::ParentPoint v_pt = atlas.parentPoint( v );
            bool at_origin = true;
            iterateGroups( pd, [&]( const size_t expanded_start, const size_t, const param::CoordinateSystem& ) {
                if( v_pt.mBaryCoordIsZero.at( expanded_start ) ) at_origin = false;
            } );
            if( at_origin ) out.emplace( v );
            return not at_origin;
        } );

        return out.value();
    }

    SmallVector<ParentPoint, 4> getFrame( const ParametricAtlas& atlas, const topology::Cell& f, const bool reverse_dart = false )
    {
        if( atlas.cmap().dim() < 2 ) throw std::runtime_error( "Cannot get face frame from atlas with dimension less than 2" );

        const param::ParentPoint ppt00 = atlas.parentPoint( f.dart() );
        const param::ParentPoint ppt10 = atlas.parentPoint( phi( atlas.cmap(), 1, f.dart() ).value() );
        const topology::Vertex third_corner(
            phi( atlas.cmap(), reverse_dart ? std::vector<int>{ 1, 1 } : std::vector<int>{ -1 }, f.dart() ).value() );
        const param::ParentPoint ppt01 = atlas.parentPoint( third_corner );

        if( atlas.cmap().dim() == 2 )
        {
            return reverse_dart ? SmallVector<ParentPoint, 4>{ { ppt10, ppt00, ppt01 } }
                                : SmallVector<ParentPoint, 4>{ { ppt00, ppt10, ppt01 } };
        }
        else
        {
            const topology::Vertex fourth_corner(
                phi( atlas.cmap(), reverse_dart ? std::vector<int>{ 2, -1 } :std::vector<int>{ 2, 1, 1 }, f.dart() ).value() );
            const param::ParentPoint ppt11 = atlas.parentPoint( fourth_corner );

            return reverse_dart ? SmallVector<ParentPoint, 4>{ { ppt10, ppt00, ppt01, ppt11 } }
                                : SmallVector<ParentPoint, 4>{ { ppt00, ppt10, ppt01, ppt11 } };
        }
    }

    SmallVector<std::pair<size_t, bool>, 3> coordinateTransform( const ParametricAtlas& atlas, const topology::Cell& f )
    {
        if( f.dim() != atlas.cmap().dim() - 1 )
        {
            throw std::runtime_error( "coordinateTransform only supports cells of dimension dim-1" );
        }
        const auto frame = getFrame( atlas, f );
        const size_t dim = atlas.cmap().dim();
        const auto maybe_phidim = phi( atlas.cmap(), dim, f.dart() );
        if( not maybe_phidim.has_value() )
        {
            throw std::runtime_error( "Face dart has no phi_dim neighbor" );
        }
        const auto other_frame = getFrame( atlas, topology::Face( maybe_phidim.value() ), true );

        const auto find_change = [dim]( const param::ParentPoint& pt, const param::ParentPoint& origin ) -> std::pair<size_t, bool> {
            for( size_t i = 0; i < dim; i++ )
            {
                if( not util::equals( pt.mPoint(i), origin.mPoint(i), 1e-10 ) )
                    return {i, copysign( 1, pt.mPoint(i) - origin.mPoint(i) ) > 0};
            }
            throw std::runtime_error( "No changing coordinate found between points" );
        };

        SmallVector<std::pair<size_t, bool>, 3> transform( dim, {0, false} );
        for( size_t coord = 0; coord < dim; coord++ )
        {
            const bool factor = ( coord != dim - 1 );
            const auto coord_1 = find_change( frame.at( coord + 1 ), frame.at( 0 ) );
            const auto coord_2 = find_change( other_frame.at( coord + 1 ), other_frame.at( 0 ) );
            transform.at( coord_1.first ) = { coord_2.first, coord_1.second == ( coord_2.second == factor ) };
        }

        return transform;
    }
}