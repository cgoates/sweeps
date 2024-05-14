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
    class ParentBasis;

    /// @brief Defines a parent basis on each element of a parametric atlas.
    ///
    class BasisComplex
    {
        public:
        /// @brief The parametric atlas associated with the basis complex.
        ///
        virtual const param::ParametricAtlas& parametricAtlas() const = 0;

        /// @brief Gives the parent basis on the given element.
        /// Only guaranteed to be implemented for element, not general cells.
        virtual const ParentBasis& parentBasis( const topology::Cell& ) const = 0;
    };
}