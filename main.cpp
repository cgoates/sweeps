#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <Laplace.hpp>
#include <cgogn/core/functions/traversals/vertex.h>
#include <VTKOutput.hpp>
#include <Logging.hpp>

cgogn::CMap3::Face nthFaceWithVertsInSet( const cgogn::CMap3& map, const std::set<VertexId>& set, const size_t n )
{
    cgogn::CMap3::Face out;
    const auto contains = [&]( const VertexId& vid ) {
        return std::find( set.begin(), set.end(), vid ) != set.end();
    };
    size_t i = 0;
    foreach_cell( map, [&]( cgogn::CMap3::Face f ) {
        out = f;
        const cgogn::Dart& d = f.dart_;
        const VertexId vid1 = index_of( map, cgogn::CMap3::Vertex( d ) );
        const VertexId vid2 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, d ) ) );
        const VertexId vid3 = index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, d ) ) );
        if( contains( vid1 ) and contains( vid2 ) and contains( vid3 ) ) i++;
        return i <= n;
    } );
    return out;
}

int main()
{
    // const SweepInput sweep_input = SweepInputTestCases::twelveTetCube();//io::loadINPFile( "/Users/caleb/sweeps/attempt-sweep/test/simple_mesh.inp", "Surface1", "Surface28" );
    const SweepInput sweep_input = io::loadINPFile( "/Users/caleb/Downloads/TorusPipe3.inp", "Surface3", "Surface4" );
    cgogn::CMap3 map;
    mapFromInput( sweep_input, map );

    const std::vector<Normal> normals = faceNormals( map );
    const Eigen::VectorXd ans = solveLaplaceSparse( map, sweep_input.zero_bcs, sweep_input.one_bcs, normals );

    const Eigen::MatrixX3d grad = gradients( map, ans, normals );
    cgogn::CMap3::Face start_face = nthFaceWithVertsInSet( map, sweep_input.zero_bcs, 10 );
    const SimplicialComplex field_line = traceField( map, start_face, centroid( map, start_face ), grad, normals );
    std::cout << "field_line: " << field_line.points << std::endl << std::endl;

    io::VTKOutputObject output( sweep_input.mesh );
    output.addVertexField( "laplace", ans );
    output.addCellField( "gradients", grad );

    // std::cout << ans.transpose() << std::endl;
    io::outputSimplicialFieldToVTK( output, "test2.vtu" );
    io::outputSimplicialFieldToVTK( map, ans, "test.vtu" );

    io::VTKOutputObject line_output( field_line );
    io::outputSimplicialFieldToVTK( line_output, "line.vtu" );

    std::cout << "dimension: " << (int)map.dimension << std::endl;
    std::cout << "Simplicial? " << cgogn::is_simplicial( map ) << std::endl;

    std::cout << "Tets: " << cgogn::nb_cells<cgogn::CMap3::Volume>( map ) << std::endl;
    std::cout << "Tris: " << cgogn::nb_cells<cgogn::CMap3::Face>( map ) << std::endl;
    std::cout << "Edges: " << cgogn::nb_cells<cgogn::CMap3::Edge>( map ) << std::endl;
    std::cout << "Vertices: " << cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) << std::endl;
}