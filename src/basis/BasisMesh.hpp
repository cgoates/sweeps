#pragma once

namespace topology
{
    class Cell;
}

namespace param
{
    class ParametricAtlas;
}

namespace basis
{
    class LocalBasis;

    /// @brief Defines a local basis on each element of a parametric atlas.
    ///
    class BasisMesh
    {
        public:
        /// @brief The parametric atlas associated with the basis mesh.
        ///
        virtual const param::ParametricAtlas& parametricAtlas() const = 0;

        /// @brief Gives the local basis on the given element.
        /// Only guaranteed to be implemented for element, not general cells.
        virtual const LocalBasis& localBasis( const topology::Cell& ) const = 0;
    };
}