#include <VTKOutput.hpp>
#include <fstream>
#include <Simplex.hpp>
#include <cassert>
#include <sstream>
#include <SplineSpace.hpp>
#include <BasisComplex.hpp>
#include <ParametricAtlas.hpp>
#include <CombinatorialMap.hpp>
#include <CombinatorialMapMethods.hpp>
#include <ParentBasis.hpp>
#include <IndexOperations.hpp>

namespace io
{
    inline std::string pieceHeader( const size_t n_points, const size_t n_cells )
    {
        std::ostringstream ss;
        ss <<R"STRING(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
  <UnstructuredGrid>
    <FieldData>
    </FieldData>
    <Piece NumberOfPoints=")STRING"
             << n_points << R"STRING(" NumberOfCells=")STRING" << n_cells << "\">";
        return ss.str();
    }

    inline std::string pieceFooter()
    {
        return "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>";
    }

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

        file << pieceHeader( n_points, n_simplices ) << "\n<CellData>";

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
)STRING" << pieceFooter();

        file.close();
    }

    SmallVector<size_t, 3> degreesOfParentBasis( const basis::ParentBasis& pb )
    {
        SmallVector<size_t, 3> out;
        for( const auto& g : pb.mBasisGroups )
        {
            for( const auto& p : g.degrees )
                out.push_back( p );
        }
        return out;
    }

    void outputBezierMeshToVTK( const basis::SplineSpace& ss,
                                const Eigen::MatrixX3d& geom,
                                const std::string& filename )
    {
        const topology::CombinatorialMap& cmap = ss.basisComplex().parametricAtlas().cmap();
        const size_t param_dim = cmap.dim();
        const size_t n_cells = cellCount( cmap, param_dim );

        std::ostringstream degrees_string;
        std::ostringstream points_string;
        std::ostringstream connectivity_string;
        std::ostringstream offsets_string;
        size_t n_points = 0;
        iterateCellsWhile( cmap, param_dim, [&]( const topology::Cell& c ) {
            SmallVector<size_t, 3> degrees = degreesOfParentBasis( ss.basisComplex().parentBasis( c ) );
            SmallVector<size_t, 3> tp_lengths;
            for( const size_t&  p : degrees ) tp_lengths.push_back( p + 1 );
            while( degrees.size() < 3 ) degrees.push_back( 0 );
            degrees_string << degrees.at( 0 ) << " " << degrees.at( 1 ) << " " << degrees.at( 2 ) << std::endl;

            const Eigen::MatrixXd ex_op = ss.extractionOperator( c );
            const std::vector<basis::FunctionId> connect = ss.connectivity( c );
            const Eigen::MatrixX3d bez_points = ex_op.transpose() * geom( connect, Eigen::all );
            
            util::iterateVTKTPOrdering( tp_lengths, [&] ( const size_t row_ii ) {
                points_string << bez_points( row_ii, 0 ) << " " << bez_points( row_ii, 1 ) << " " << bez_points( row_ii, 2 ) << std::endl;
            } );

            for( Eigen::Index i = 0; i < bez_points.rows(); i++ )
                connectivity_string << ( n_points + i ) << " ";
            n_points += bez_points.rows();
            offsets_string << n_points << " ";

            return true;
        } );

        std::ofstream file;
        file.open( filename );

        file << pieceHeader( n_points, n_cells ) << "\n<CellData HigherOrderDegrees=\"degrees\">";

        file << R"STRING(
        <DataArray type="UInt8" Name="degrees" NumberOfComponents="3" format="ascii">)STRING" << std::endl;

        file << degrees_string.str();

        file << R"STRING(
        </DataArray>
      </CellData>
      <PointData>
      </PointData>
      <Points>
        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">
)STRING";

        file << points_string.str();

        file << R"STRING(
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">
)STRING";

        file << connectivity_string.str();

        file << R"STRING(
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
)STRING";

        file << offsets_string.str();

        file << R"STRING(
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)STRING";

        iterateCellsWhile( cmap, param_dim, [&]( const auto& ) {
            switch( param_dim )
            {
                // See documentation at
                // https://www.kitware.com/main/wp-content/uploads/2020/03/Implementation-of-rational-Be%CC%81zier-cells-into-VTK-Report.pdf
                // page 3
                case 0: file << 1 << " "; break;
                case 1: file << 75 << " "; break;
                case 2: file << 77 << " "; break;
                case 3: file << 79 << " "; break;
                default: break;
            }
            return true;
        } );

        file << R"STRING(
        </DataArray>
      </Cells>
)STRING" << pieceFooter();

        file.close();
    }
} // namespace io