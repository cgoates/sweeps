#include <ParametricAtlas.hpp>
#include <CombinatorialMapMethods.hpp>

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
            iterateGroups( pd, [&]( const size_t expanded_start, const size_t explicit_start, const param::CoordinateSystem& cs ) {
                if( v_pt.mBaryCoordIsZero.at( expanded_start ) ) at_origin = false;
            } );
            if( at_origin ) out.emplace( v );
            return not at_origin;
        } );

        return out.value();
    }
}