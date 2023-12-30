#include <VTKOutput.hpp>
#include <fstream>
#include <Simplex.hpp>
#include <cgogn/core/types/cell_marker.h>
#include <cgogn/core/types/maps/cmap/cmap3.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/io/volume/volume_import.h>

namespace io
{

Simplex volumeConnectivity( const cgogn::CMap3& map, const cgogn::CMap3::Volume& v )
{
    return Simplex( VertexId( cgogn::index_of( map, cgogn::CMap3::Vertex( v.dart_ ) ) ),
                    VertexId( cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi1( map, v.dart_ ) ) ) ),
                    VertexId( cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi_1( map, v.dart_ ) ) ) ),
                    VertexId( cgogn::index_of( map, cgogn::CMap3::Vertex( cgogn::phi<2,-1>( map, v.dart_ ) ) ) ) );
}

void outputSimplicialFieldToVTK( const cgogn::CMap3& map, const Eigen::MatrixXd& data, const std::string& filename )
{
    assert( cgogn::is_simplicial( map ) );
    assert( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ) == data.rows() );

    const size_t n_simplices = cgogn::nb_cells<cgogn::CMap3::Volume>( map );

    std::ofstream file;
    file.open( filename );

    file << R"STRING(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
  <UnstructuredGrid>
    <FieldData>
    </FieldData>
    <Piece NumberOfPoints=")STRING" << data.rows()
    << R"STRING(" NumberOfCells=")STRING"
    << n_simplices << R"STRING(">
      <CellData>
      </CellData>
      <PointData>
        <DataArray type="Float64" Name="DISPL" NumberOfComponents=")STRING" 
        << data.cols() << R"STRING(" format="ascii">)STRING" << std::endl;

    for( Eigen::Index row = 0; row < data.rows(); row++ )
    {
        for( Eigen::Index col = 0; col < data.cols(); col++ )
        {
            file << data( row, col ) << " ";
        }
    }

    file << R"STRING(
        </DataArray>
      </PointData>
      <Points>
        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">
)STRING";

    const auto position = cgogn::get_attribute<Eigen::Vector3d, cgogn::CMap3::Vertex>( map, "position" );
    Eigen::MatrixX3d points( cgogn::nb_cells<cgogn::CMap3::Vertex>( map ), 3 );

    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Vertex v ){
        const auto vid = cgogn::index_of( map, v );
        points.row( vid ) = cgogn::value<Eigen::Vector3d>( map, position, v );
        return true;
    } );

    for( Eigen::Index row = 0; row < points.rows(); row++ )
    {
        for( Eigen::Index col = 0; col < points.cols(); col++ )
        {
            file << points( row, col ) << " ";
        }
    }

    file << R"STRING(
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
)STRING";

    cgogn::foreach_cell( map, [&]( cgogn::CMap3::Volume v ) {
        const Simplex simplex = volumeConnectivity( map, v );
        for( size_t vert_ii = 0; vert_ii <= simplex.dim(); vert_ii++ )
        {
            file << simplex.vertex( vert_ii ).id() << " ";
        }
        return true;
    } );

    file << R"STRING(
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
)STRING";

    for( size_t simplex_ii = 0; simplex_ii < n_simplices; simplex_ii++ )
    {
        file << ( ( simplex_ii + 1 ) * 4 ) << " ";
    }
    file << R"STRING(
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)STRING";

    for( size_t simplex_ii = 0; simplex_ii < n_simplices; simplex_ii++ )
    {
        file << 10 << " ";
    }
    
    file << R"STRING(
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)STRING";

    file.close();
}

}