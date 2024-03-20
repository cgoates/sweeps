#pragma once
#include <string>
#include <Eigen/Dense>
#include <SweepInput.hpp>
#include <map>

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
} // namespace io