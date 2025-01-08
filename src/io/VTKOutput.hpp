#pragma once
#include <string>
#include <CustomEigen.hpp>
#include <SimplicialComplex.hpp>
#include <VertexPositionsFunc.hpp>
#include <map>

namespace topology
{
    class CombinatorialMap;
    class Edge;
    class Cell;
}

namespace basis
{
    class SplineSpace;
};

namespace io
{
    class VTKOutputObject
    {
        public:
        VTKOutputObject( const SimplicialComplex& complex ) : mComplex( complex ) {}

        void addVertexField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field );
        void addCellField( const std::string& label, const Eigen::Ref<const Eigen::MatrixXd> field );

        const SimplicialComplex& complex() const { return mComplex; }
        const std::map<std::string, Eigen::MatrixXd>& vertexFields() const { return mVertexFields; }
        const std::map<std::string, Eigen::MatrixXd>& cellFields() const { return mCellFields; }

        private:
        const SimplicialComplex& mComplex;
        std::map<std::string, Eigen::MatrixXd> mVertexFields;
        std::map<std::string, Eigen::MatrixXd> mCellFields;
    };

    void outputSimplicialFieldToVTK( const VTKOutputObject& output, const std::string& filename );

    /// @brief write out a spline as bezier cells to vtk
    /// @param ss  The spline space
    /// @param geom The control points, in a n_ctrl_pts x spatial_dim sized matrix
    /// @param filename The output filename, including the extension.
    void outputBezierMeshToVTK( const basis::SplineSpace& ss,
                                const MatrixX3dMax& geom,
                                const std::string& filename );

    void outputPartialBezierMeshToVTK( const basis::SplineSpace& ss,
                                       const MatrixX3dMax& geom,
                                       const std::string& filename,
                                       const std::function<void( const std::function<void( const topology::Cell& )>& )>& cell_iterator );

    void outputEdges( const topology::CombinatorialMap& cmap,
                      const VertexPositionsFunc& positions,
                      const std::vector<topology::Edge>& edges,
                      const std::string& filename );

    void outputCMap( const topology::CombinatorialMap& cmap,
                     const VertexPositionsFunc& positions,
                     const std::string& filename );

    void outputDualFace( const topology::CombinatorialMap& cmap,
                         const VertexPositionsFunc& positions,
                         const topology::Edge& e,
                         const std::string& postfix );
} // namespace io