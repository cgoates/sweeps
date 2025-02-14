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
        return 0;// Added to make linux build happy
    }

    size_t numFunctions( const ParentBasis& pb )
    {
        return std::accumulate( pb.mBasisGroups.begin(), pb.mBasisGroups.end(), 1, [&]( const size_t accum, const BarycentricBasis& bb ) {
            return accum * numFunctions( bb );
        } );
    }

    SmallVector<size_t, 3> degrees( const ParentBasis& pb )
    {
        SmallVector<size_t, 3> degrees;
        for( const auto& group : pb.mBasisGroups )
        {
            for( const auto& p : group.degrees )
                degrees.push_back( p );
        }
        return degrees;
    }

    size_t numVectorComponents( const ParentBasis& pb )
    {
        if( std::any_of( pb.mBasisGroups.begin(), pb.mBasisGroups.end(), []( const BarycentricBasis& bb ) {
                return bb.type == BasisType::DivConformingBernstein;
            } ) )
        {
            if( pb.mBasisGroups.size() == 1 )
                return dim( pb.mParentDomain );
            else
                throw std::runtime_error( "Cannot tensor product vector valued bases" );
        }
        return 1;
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

    std::ostream& operator<<( std::ostream& o, const BarycentricBasis& bb )
    {
        o << "BarycentricBasis( ";
        switch( bb.type )
        {
            case BasisType::Bernstein:
                o << "Bernstein";
                break;
            case BasisType::DivConformingBernstein:
                o << "DivConformingBernstein";
        }
        o << ", " << bb.degrees << " )";
        return o;
    }

    std::ostream& operator<<( std::ostream& o, const ParentBasis& pb )
    {
        o << "ParentBasis( " << pb.mParentDomain << ", " << pb.mBasisGroups << " )";
        return o;
    }

    SmallVector<ParentBasis, 3> componentBases( const ParentBasis& pb )
    {
        if( numVectorComponents( pb ) == 1 ) return { pb };

        SmallVector<ParentBasis, 3> out;
        for( const auto& group : pb.mBasisGroups )
        {
            if( group.type == BasisType::DivConformingBernstein )
            {
                SmallVector<BarycentricBasis, 3> primal_groups;
                std::transform( group.degrees.begin(), group.degrees.end(), std::back_inserter( primal_groups ), []( const size_t degree ) {
                    return BarycentricBasis{ BasisType::Bernstein, { degree - 1 } };
                } );
                for( size_t i = 0; i < group.degrees.size(); i++ )
                {
                    SmallVector<BarycentricBasis, 3> primal_groups_i = primal_groups;
                    primal_groups_i.at( i ) = BarycentricBasis{ BasisType::Bernstein, { group.degrees.at( i ) } };
                    out.push_back( ParentBasis{ pb.mParentDomain, primal_groups_i } );
                }
            }
            else
            {
                throw std::runtime_error( "Cannot tensor product div conforming and scalar bases" );
            }
        }

        return out;
    }
}