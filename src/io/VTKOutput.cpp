#include <VTKOutput.hpp>
#include <fstream>
#include <Simplex.hpp>
#include <cassert>

namespace io
{
    void VTKOutputObject::addVertexField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field )
    {
        mVertexFields.emplace( label, field );
    }

    void VTKOutputObject::addCellField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field )
    {
        mCellFields.emplace( label, field );
    }

    void outputSimplicialFieldToVTK( const VTKOutputObject& to_output, const std::string& filename )
    {
        const size_t n_simplices = to_output.complex().simplices.size();
        const size_t n_points = to_output.complex().points.size();

        std::ofstream file;
        file.open( filename );

        file << R"STRING(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
  <UnstructuredGrid>
    <FieldData>
    </FieldData>
    <Piece NumberOfPoints=")STRING"
             << n_points << R"STRING(" NumberOfCells=")STRING" << n_simplices << R"STRING(">
      <CellData>)STRING";

        for( const auto& [label, data] : to_output.cellFields() )
        {
            file << R"STRING(
        <DataArray type="Float64" Name=")STRING" << label << R"STRING(" NumberOfComponents=")STRING"
                 << data.cols() << R"STRING(" format="ascii">)STRING" << std::endl;

            for( Eigen::Index row = 0; row < data.rows(); row++ )
            {
                for( Eigen::Index col = 0; col < data.cols(); col++ )
                {
                    file << data( row, col ) << " ";
                }
            }

            file << R"STRING(
        </DataArray>)STRING";
        }

        file << R"STRING(
      </CellData>
      <PointData>)STRING";

        for( const auto& [label, data] : to_output.vertexFields() )
        {
            file << R"STRING(
        <DataArray type="Float64" Name=")STRING" << label << R"STRING(" NumberOfComponents=")STRING"
                 << data.cols() << R"STRING(" format="ascii">)STRING" << std::endl;

            for( Eigen::Index row = 0; row < data.rows(); row++ )
            {
                for( Eigen::Index col = 0; col < data.cols(); col++ )
                {
                    file << data( row, col ) << " ";
                }
            }

            file << R"STRING(
        </DataArray>)STRING";
        }

        file << R"STRING(
      </PointData>
      <Points>
        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">
)STRING";

        for( const Eigen::Vector3d& point : to_output.complex().points )
        {
            file << point( 0 ) << " " << point( 1 ) << " " << point( 2 ) << " ";
        }

        file << R"STRING(
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
)STRING";

        for( const Simplex& simplex : to_output.complex().simplices )
        {
            for( size_t vert_ii = 0; vert_ii <= simplex.dim(); vert_ii++ )
            {
                file << simplex.vertex( vert_ii ).id() << " ";
            }
        }

        file << R"STRING(
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
)STRING";

        size_t offset = 0;
        for( const Simplex& simplex : to_output.complex().simplices )
        {
            offset += simplex.dim() + 1;
            file << offset << " ";
        }
        file << R"STRING(
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)STRING";

        for( const Simplex& simplex : to_output.complex().simplices )
        {
            switch( simplex.dim() )
            {
                // See documentation at https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf page 9
                case 0: file << 1 << " "; break;
                case 1: file << 3 << " "; break;
                case 2: file << 5 << " "; break;
                case 3: file << 10 << " "; break;
                default: break;
            }
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

} // namespace io