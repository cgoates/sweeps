#include <LocalBasis.hpp>
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

    LocalBasis bernsteinSimplex( const size_t dim, const size_t degree )
    {
        return LocalBasis{ param::simplexDomain( dim ), { bernsteinBasis( degree ) } };
    }
    
    LocalBasis bernsteinCube( const size_t dim, const size_t degree )
    {
        return LocalBasis{ param::cubeDomain( dim ), { dim, bernsteinBasis( degree ) } };
    }
}