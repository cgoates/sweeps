#include <CoonsPatch.hpp>
#include <IndexOperations.hpp>

namespace fitting
{
    Eigen::VectorXd nLinearInterpolation( const Eigen::VectorXd& pt, const std::vector<Eigen::VectorXd>& data )
    {
        if( pow( 2, pt.size() ) != data.size() )
            throw std::runtime_error( "Invalid inputs to nLinearInterpolation");

        std::vector<Eigen::VectorXd> working_copy = data;

        size_t offset = data.size() / 2;
        for( int i = pt.size() - 1; i >= 0; i-- )
        {
            for( size_t j = 0; j < offset; j++ )
            {
                working_copy.at( j ) =
                    ( ( 1 - pt( i ) ) * working_copy.at( j ) + pt( i ) * working_copy.at( j + offset ) ).eval();
            }
            offset /= 2;
        }
        return working_copy.front();
    }

    Eigen::VectorXd coonsInterpolation( const util::IndexVec& tp_id,
                                        const SmallVector<Eigen::VectorXd, 3>& linear_coeffs,
                                        const util::IndexVec& num_funcs,
                                        const size_t data_dim,
                                        const std::vector<Eigen::MatrixXd>& boundary_cpts )
    {
        const size_t dim = linear_coeffs.size();

        const Vector3dMax point = [&linear_coeffs, &tp_id]() {
            Vector3dMax out( tp_id.size() );
            for( size_t i = 0; i < tp_id.size(); i++ ) out( i ) = linear_coeffs.at( i )( tp_id.at( i ) );
            return out;
        }();

        constexpr bool to_zero = false;
        const auto clamped_cpt = [&]( const util::IndexVec& tp_id, const SmallVector<size_t, 3>& index, const SmallVector<bool, 3>& clamp_positive ) {
            util::IndexVec clamped_id = tp_id;
            for( size_t i = 0; i < index.size() and i < clamp_positive.size(); i++ )
            {
                clamped_id.at( index.at( i ) ) = clamp_positive.at( i ) ? num_funcs.at( index.at( i ) ) - 1 : 0;
            }

            Eigen::VectorXd out = Eigen::VectorXd::Zero( data_dim );
            // If we are dropping multiple indices, there are multiple boundary coefficient sets that we could grab the
            // clamped control point from, so here we just average these points. Typically these will all be the same,
            // but because we are averaging this is not required; lower dimensional boundaries are just defined as the
            // average of the provided values.
            for( size_t i = 0; i < index.size() and i < clamp_positive.size(); i++ )
            {
                const size_t which_boundary = 2 * index.at( i ) + ( clamp_positive.at( i ) ? 1 : 0 );
                out += boundary_cpts.at( which_boundary )
                           .row( util::flatten( util::dropIndex( clamped_id, index.at( i ) ),
                                                util::dropIndex( num_funcs, index.at( i ) ) ) ) /
                       index.size();
            }

            return out;
        };

        Eigen::VectorXd result = Eigen::VectorXd::Zero( data_dim );

        // Linear interpolations
        for( size_t dim_ii = 0; dim_ii < dim; dim_ii++ )
        {
            result += nLinearInterpolation(
                Vector1d( point( dim_ii ) ),
                { clamped_cpt( tp_id, { dim_ii }, { to_zero } ), clamped_cpt( tp_id, { dim_ii }, { not to_zero } ) } );
        }

        // Bilinear interpolations
        for( size_t dim_ii = 0; dim_ii < dim - 1; dim_ii++ )
        {
            for( size_t dim_jj = dim_ii + 1; dim_jj < dim; dim_jj++ )
            {
                result -= nLinearInterpolation(
                    Eigen::Vector2d( point( dim_ii ), point( dim_jj ) ),
                    { clamped_cpt( tp_id, { dim_ii, dim_jj }, { to_zero,     to_zero     } ),
                      clamped_cpt( tp_id, { dim_ii, dim_jj }, { not to_zero, to_zero     } ),
                      clamped_cpt( tp_id, { dim_ii, dim_jj }, { to_zero,     not to_zero } ),
                      clamped_cpt( tp_id, { dim_ii, dim_jj }, { not to_zero, not to_zero } ) } );
            }
        }

        if( dim > 2 )
        {
            // Trilinear interpolation
            result += nLinearInterpolation(
                point,
                { clamped_cpt( tp_id, { 0, 1, 2 }, { to_zero,     to_zero,     to_zero     } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { not to_zero, to_zero,     to_zero     } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { to_zero,     not to_zero, to_zero     } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { not to_zero, not to_zero, to_zero     } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { to_zero,     to_zero,     not to_zero } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { not to_zero, to_zero,     not to_zero } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { to_zero,     not to_zero, not to_zero } ),
                  clamped_cpt( tp_id, { 0, 1, 2 }, { not to_zero, not to_zero, not to_zero } ) } );
        }

        return result;
    }

    Eigen::MatrixXd coonsPatch( const basis::TPSplineSpace& ss, const std::vector<Eigen::MatrixXd>& boundary_cpts )
    {
        const std::vector<std::reference_wrapper<const basis::BSplineSpace1d>> spline_1ds =
            constituentSplines( ss ).value_or( std::vector<std::reference_wrapper<const basis::BSplineSpace1d>>() );
        
        if( spline_1ds.size() == 0 ) throw std::runtime_error( "Cannot create a coons patch on TPSplineSpace with non-1d constituents." );
        if( boundary_cpts.size() != 2 * spline_1ds.size() ) throw std::runtime_error( "Incorrect number of boundary coefficient sets specified." );

        const size_t data_dim = boundary_cpts.at( 0 ).cols();
        SmallVector<Eigen::VectorXd, 3> linear_coeffs;
        util::IndexVec num_funcs;

        for( const basis::BSplineSpace1d& spline : spline_1ds )
        {
            linear_coeffs.push_back(
                grevillePoints( spline.knotVector(), degrees( spline.basisComplex().defaultParentBasis() ).at( 0 ) ) /
                spline.knotVector().knot( spline.knotVector().size() - 1 ) );
            num_funcs.push_back( spline.numFunctions() );
        }

        Eigen::MatrixXd coons_coefficients = Eigen::MatrixXd::Zero( ss.numFunctions(), data_dim );
        util::iterateTensorProduct( num_funcs, [&]( const util::IndexVec& tp_id ) {
            const size_t func_id = util::flatten( tp_id, num_funcs );
            coons_coefficients.row( func_id ) = coonsInterpolation( tp_id, linear_coeffs, num_funcs, data_dim, boundary_cpts );
        } );

        return coons_coefficients;
    }
}