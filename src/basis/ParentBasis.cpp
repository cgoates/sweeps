#include <ParentBasis.hpp>
#include <ParentDomain.hpp>

namespace basis
{
    size_t degree( const BarycentricBasis& bb )
    {
        return bb.num_functions - 1;
    }

    BarycentricBasis bernsteinBasis( const size_t degree )
    {
        return BarycentricBasis{ BasisType::Bernstein, degree + 1 };
    }

    ParentBasis bernsteinSimplex( const size_t dim, const size_t degree )
    {
        return ParentBasis{ param::simplexDomain( dim ), { bernsteinBasis( degree ) } };
    }
    
    ParentBasis bernsteinCube( const size_t dim, const size_t degree )
    {
        return ParentBasis{ param::cubeDomain( dim ), { dim, bernsteinBasis( degree ) } };
    }
}