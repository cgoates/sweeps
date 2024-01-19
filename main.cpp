#include <SimplexUtilities.hpp>
#include <AbaqusInput.hpp>
#include <Laplace.hpp>
#include <cgogn/core/functions/traversals/vertex.h>
#include <VTKOutput.hpp>
#include <Logging.hpp>

void foreachFaceWithVertsInSet( const cgogn::CMap3& map,
                                const std::set<VertexId>& set,
                                const std::function<bool( const cgogn::CMap3::Face&, const size_t n )>& callback )
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
        if( contains( vid1 ) and contains( vid2 ) and contains( vid3 ) ) return callback( f, i++ );
        return true;
    } );
}


void append( SimplicialComplex& append_to, const SimplicialComplex& to_append )
{
    const size_t offset = append_to.points.size();
    append_to.points.insert( append_to.points.end(), to_append.points.begin(), to_append.points.end() );
    std::transform( to_append.simplices.cbegin(), to_append.simplices.cend(), std::back_inserter( append_to.simplices ), [&]( const Simplex& s ) {
        // FIXME: NOT GENERAL!!!
        return Simplex( s.vertex( 0 ).id() + offset, s.vertex( 1 ).id() + offset );
    } );
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

    io::VTKOutputObject output( sweep_input.mesh );
    output.addVertexField( "laplace", ans );
    output.addCellField( "gradients", grad );

    // std::cout << ans.transpose() << std::endl;
    io::outputSimplicialFieldToVTK( output, "test.vtu" );

    SimplicialComplex all_lines;

    foreachFaceWithVertsInSet( map, sweep_input.zero_bcs, [&]( const cgogn::CMap3::Face& start_face, const size_t i ){
        const SimplicialComplex field_line = traceField( map, start_face, centroid( map, start_face ), grad, normals );
        // std::cout << "field_line: " << field_line.points << std::endl << std::endl;
        if( field_line.simplices.size() < 40 )std::cout << "Line: " << i << " Length: " << field_line.simplices.size() << std::endl;
        append( all_lines, field_line );
        return true;
    } );

    io::VTKOutputObject line_output( all_lines );
    io::outputSimplicialFieldToVTK( line_output, "line_full.vtu" );

    foreachFaceWithVertsInSet( map, sweep_input.zero_bcs, [&]( const cgogn::CMap3::Face& start_face, const size_t i ){
        if( i == 184 )
        {
            const SimplicialComplex field_line = traceField( map, start_face, centroid( map, start_face ), grad, normals, true );
            // std::cout << "field_line: " << field_line.points << std::endl << std::endl;
        }
        return true;
    } );

    std::cout << "dimension: " << (int)map.dimension << std::endl;
    std::cout << "Simplicial? " << cgogn::is_simplicial( map ) << std::endl;

    std::cout << "Tets: " << cgogn::nb_cells<cgogn::CMap3::Volume>( map ) << std::endl;
    std::cout << "Tris: " << cgogn::nb_cells<cgogn::CMap3::Face>( map ) << std::endl;
    std::cout << "Edges: " << cgogn::nb_cells<cgogn::CMap3::Edge>( map ) << std::endl;
    std::cout << "Vertices: " << cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) << std::endl;
}