#include<SimplexUtilities.hpp>
#include<AbaqusInput.hpp>
#include<Laplace.hpp>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include<VTKOutput.hpp>

int main()
{
    // const SweepInput sweep_input = twoTets();
    // const SweepInput sweep_input = twelveTetCube();
    // const SweepInput sweep_input = refinedCube( {3,4, 2} );
    // const SweepInput sweep_input = bentRefinedCube( {3, 3, 30 } );
    // const SweepInput sweep_input = loadINPFile( "/Users/caleb/sweeps/attempt-sweep/test/simple_mesh.inp", "Surface1", "Surface28" );
    const SweepInput sweep_input = io::loadINPFile( "/Users/caleb/Downloads/TorusPipe.inp", "Surface3", "Surface4" );
    cgogn::CMap3 map;
    SimplexUtilities::mapFromInput( sweep_input, map );

    /*
        TODO:
      X - build laplace operator
      X - solve laplace equation for harmonic function
      X - Tests of cgogn import
        - Tests of Laplace solve
        - Mesh input from gmsh with floor/ceiling specs
    */

    // const Eigen::VectorXd ans1 = solveLaplace( map, sweep_input.zero_bcs, sweep_input.one_bcs );
    const Eigen::VectorXd ans = solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs );

    std::cout << ans.transpose() << std::endl;
    io::outputSimplicialFieldToVTK( map, ans, "test.vtu" );

    std::cout << "dimension: " << (int)map.dimension << std::endl;
    std::cout << "Simplicial? " << cgogn::is_simplicial( map ) << std::endl;

    std::cout << "Tets: " << cgogn::nb_cells<cgogn::CMap3::Volume>( map ) << std::endl;
    std::cout << "Tris: " << cgogn::nb_cells<cgogn::CMap3::Face>( map ) << std::endl;
    std::cout << "Edges: " << cgogn::nb_cells<cgogn::CMap3::Edge>( map ) << std::endl;
    std::cout << "Vertices: " << cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) << std::endl;
}