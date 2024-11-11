#include <SweepAPIMethods.hpp>
#include <SweepInput.hpp>
#include <Foliation.hpp>
#include <VTKOutput.hpp>
#include <TriangleMeshMapping.hpp>
#include <SimplexUtilities.hpp>
#include <CombinatorialMapMethods.hpp>

namespace api
{
    void outputLevelSetsAndTraces( const Sweep& sweep,
                                   const std::vector<double>& level_set_values,
                                   const std::vector<Eigen::Vector2d>& trace_points,
                                   const std::string& output_prefix )
    {
        if( sweep.source.size() == 0 or sweep.target.size() == 0 )
            throw std::invalid_argument( "No source or target surface specified" );
        if( sweep.mesh.points.size() == 0 or sweep.mesh.simplices.size() == 0 )
            throw std::invalid_argument( "No tet mesh supplied" );
        if( level_set_values.size() < 2 )
            throw std::invalid_argument( "Insufficient level sets specified; please include at least two values." );

        const SweepInput sweep_input = [&sweep]() {
            const auto& m = sweep.mesh;
            std::vector<bool> zero_bcs( m.points.size(), false );
            std::vector<bool> one_bcs( m.points.size(), false );
            for( const VertexId::Type& vid : sweep.source ) zero_bcs.at( vid ) = true;
            for( const VertexId::Type& vid : sweep.target ) one_bcs.at( vid ) = true;
            return SweepInput( m, zero_bcs, one_bcs );
        }();

        reparam::levelSetBasedTracing( sweep_input, level_set_values, [&]( const std::vector<reparam::FoliationLeaf>& leaves ) {
            SimplicialComplex level_sets;
            SimplicialComplex param_out;

            for( const auto& leaf : leaves )
            {
                iterateCellsWhile( leaf.space_mapping->parametricAtlas().cmap(), 2, [&]( const topology::Face& f ) {
                    addTriangleNoDuplicateChecking( level_sets,
                                                    triangleOfFace<3>( leaf.space_mapping->parametricAtlas().cmap(),
                                                                    leaf.space_mapping->vertPositions(),
                                                                    f ) );
                    return true;
                } );
            }

            // Draw some lines!
            param_out.points.reserve( trace_points.size() * level_set_values.size() );
            param_out.simplices.reserve( trace_points.size() * level_set_values.size() );

            for( const auto& circle_pt : trace_points )
            {
                const size_t offset = param_out.points.size();
                size_t i = 0;
                for( const auto& leaf : leaves )
                {
                    const auto& param_pt = leaf.tutte_mapping->maybeInverse( circle_pt );
                    if( not param_pt.has_value() )
                    {
                        std::cerr << "NO VALUE pt: " << circle_pt.transpose() << " level: " << i << std::endl;
                        break;
                    }
                    const auto space_pt = leaf.space_mapping->evaluate( param_pt.value().first, param_pt.value().second );
                    param_out.points.push_back( space_pt );
                    i++;
                }

                if( i == 0 ) continue;
                for( size_t simplex_ii = 0; simplex_ii < i - 1; simplex_ii++ )
                {
                    param_out.simplices.push_back( Simplex( offset + simplex_ii, offset + simplex_ii + 1 ) );
                }
            }

            io::VTKOutputObject output( level_sets );
            io::outputSimplicialFieldToVTK( output, output_prefix + "_level_sets.vtu" );

            io::VTKOutputObject output2( param_out );
            io::outputSimplicialFieldToVTK( output2, output_prefix + "_foliation_param.vtu" );
        } );
    }
}