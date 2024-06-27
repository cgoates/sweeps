#include <ParentBasis.hpp>
#include <ParentDomain.hpp>

namespace basis
{
    size_t numFunctions( const BarycentricBasis& bb )
    {
        switch( bb.type )
        {
            case BasisType::Bernstein: return bb.degrees.at( 0 ) + 1;
            case BasisType::DivConformingBernstein:
            {
                size_t result = 0;
                for( size_t i = 0; i < bb.degrees.size(); i++ )
                {
                    size_t result_i = 1;
                    for( size_t j = 0; j < bb.degrees.size(); j++ )
                        result_i *= j == i ? bb.degrees.at( j ) + 1 : bb.degrees.at( j );
                    result += result_i;
                }
                return result;
            }
        }
    }

    BarycentricBasis bernsteinBasis( const size_t degree )
    {
        return BarycentricBasis{ BasisType::Bernstein, { degree } };
    }

    BarycentricBasis divConformingBernsteinBasis( const size_t dim, const size_t degree )
    {
        return divConformingBernsteinBasis( SmallVector<size_t, 3>( dim, degree ) );
    }

    BarycentricBasis divConformingBernsteinBasis( const SmallVector<size_t, 3>& degrees )
    {
        return BarycentricBasis{ BasisType::DivConformingBernstein, degrees };
    }

    ParentBasis bernsteinSimplex( const size_t dim, const size_t degree )
    {
        return ParentBasis{ param::simplexDomain( dim ), { bernsteinBasis( degree ) } };
    }
    
    ParentBasis bernsteinCube( const size_t dim, const size_t degree )
    {
        return ParentBasis{ param::cubeDomain( dim ), { dim, bernsteinBasis( degree ) } };
    }

    ParentBasis divConformingBernsteinCube( const size_t dim, const size_t primal_degree )
    {
        return ParentBasis{ param::cubeDomain( dim ), { divConformingBernsteinBasis( dim, primal_degree ) } };
    }

    ParentBasis tensorProduct( const ParentBasis& pb1, const ParentBasis& pb2 )
    {
        SmallVector<BarycentricBasis, 3> combined_groups = pb1.mBasisGroups;
        for( const auto& group : pb2.mBasisGroups ) combined_groups.push_back( group );
        return ParentBasis{ tensorProduct( pb1.mParentDomain, pb2.mParentDomain ), combined_groups };
    }
}