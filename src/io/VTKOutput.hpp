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
    class Dart;
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


    class BezierOutputObject
    {
        public:
        BezierOutputObject( const basis::SplineSpace& ss, const MatrixX3dMax& geom ) : mSS( ss ), mGeom( geom ) {}

        const basis::SplineSpace& ss() const { return mSS; }
        const MatrixX3dMax& geom() const { return mGeom; }

        void addBezierField( const std::string& label, const basis::SplineSpace& ss, const MatrixX3dMax& geom )
        {
            mBezierFields.emplace( label, BezierField{ std::cref( ss ), geom } );
        }

        void addBezierField( const std::string& label, const MatrixX3dMax& geom )
        {
            mBezierFields.emplace( label, BezierField{ std::nullopt, geom } );
        }

        struct BezierField
        {
            std::optional<std::reference_wrapper<const basis::SplineSpace>> ss;
            MatrixX3dMax geom;
        };

        const std::map<std::string, BezierField>& bezierFields() const { return mBezierFields; }

        private:

        const basis::SplineSpace& mSS;
        const MatrixX3dMax& mGeom;
        std::map<std::string, BezierField> mBezierFields;
    };

    void outputSimplicialFieldToVTK( const VTKOutputObject& output, const std::string& filename );

    /// @brief write out a spline as bezier cells to vtk
    /// @param ss  The spline space
    /// @param geom The control points, in a n_ctrl_pts x spatial_dim sized matrix
    /// @param filename The output filename, including the extension.
    void outputBezierMeshToVTK( const basis::SplineSpace& ss,
                                const MatrixX3dMax& geom,
                                const std::string& filename );

    void outputBezierMeshToVTK( const BezierOutputObject& output, const std::string& filename );

    void outputPartialBezierMeshToVTK( const basis::SplineSpace& ss,
                                       const MatrixX3dMax& geom,
                                       const std::string& filename,
                                       const std::function<void( const std::function<void( const topology::Cell& )>& )>& cell_iterator );

    void outputPartialBezierMeshToVTK( const BezierOutputObject& output,
                                       const std::string& filename,
                                       const std::function<void( const std::function<void( const topology::Cell& )>& )>& cell_iterator );

    void outputEdgeChain( const topology::CombinatorialMap& cmap,
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

    void outputDarts( const topology::CombinatorialMap& cmap,
                      const VertexPositionsFunc& positions,
                      const std::string& filename,
                      const std::map<std::string, std::function<Eigen::VectorXd( const topology::Dart& )>>& fields = {},
                      const std::function<bool( const topology::Dart& )>& filter = []( const auto& ) { return true; } );
} // namespace io