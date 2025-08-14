#include <ParentBasisEval.hpp>
#include <ParentBasis.hpp>
#include <ParentPoint.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <Logging.hpp>

namespace eval
{
    double binomial( const int n, const int k )
    {
        double out = 1.0;
        for( int i = 1, stop = std::min( k, n - k ); i < stop; i++ ) out *= double( n + 1 - i ) / i;

        return out;
    }

    Eigen::VectorXd binomial( const size_t n )
    {
        Eigen::VectorXd out = Eigen::VectorXd::Ones( n + 1 );
        double accum = 1.0;
        for( size_t i = 1, e = n / 2; i <= e; i++ )
        {
            accum *= double( n + 1 - i ) / i;
            out( i ) = accum;
            out( n - i ) = accum;
        }
        return out;
    }

    Eigen::VectorXd bernstein( const size_t p, const double s )
    {
        const Eigen::VectorXd powers = Eigen::VectorXd::LinSpaced( p + 1, 0, p );
        const Eigen::VectorXd binom = binomial( p );
        return binom.cwiseProduct( Eigen::VectorXd(
            pow( s, powers.array() )
                .cwiseProduct( pow( 1.0 - s, ( Eigen::VectorXd::Constant( p + 1, p ) - powers ).array() ) ) ) );
    }

    Eigen::VectorXd bernsteinFirstDeriv( const size_t p, const double s )
    {
        if( p < 1 ) return Eigen::VectorXd::Zero( p + 1 );
        const Eigen::VectorXd lower_degree = bernstein( p - 1, s );
        return p * ( ( Eigen::VectorXd( p + 1 ) << 0, lower_degree ).finished() -
                     ( Eigen::VectorXd( p + 1 ) << lower_degree, 0 ).finished() );
    }

    Eigen::VectorXd bernsteinSecondDeriv( const size_t p, const double s )
    {
        if( p < 2 ) return Eigen::VectorXd::Zero( p + 1 );
        const Eigen::VectorXd lower_degree = bernstein( p - 2, s );
        return p * ( p - 1 ) *
               ( ( Eigen::VectorXd( p + 1 ) << 0, 0, lower_degree ).finished() -
                 2 * ( Eigen::VectorXd( p + 1 ) << 0, lower_degree, 0 ).finished() +
                 ( Eigen::VectorXd( p + 1 ) << lower_degree, 0, 0 ).finished() );
    }

    Eigen::VectorXd bernsteinTP( const size_t p, const size_t q, const double s, const double t )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );

        return Eigen::kroneckerProduct( t_evals, s_evals );
    }

    Eigen::MatrixX2d bernsteinTPFirstDeriv( const size_t p, const size_t q, const double s, const double t )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );
        const Eigen::VectorXd s_deriv = bernsteinFirstDeriv( p, s );
        const Eigen::VectorXd t_deriv = bernsteinFirstDeriv( q, t );

        return ( Eigen::MatrixX2d( s_evals.size() * t_evals.size(), 2 ) << Eigen::kroneckerProduct( t_evals, s_deriv ),
                 Eigen::kroneckerProduct( t_deriv, s_evals ) )
            .finished();
    }

    Eigen::MatrixX3d bernsteinTPSecondDeriv( const size_t p, const size_t q, const double s, const double t )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );
        const Eigen::VectorXd s_deriv = bernsteinFirstDeriv( p, s );
        const Eigen::VectorXd t_deriv = bernsteinFirstDeriv( q, t );
        const Eigen::VectorXd s_dderiv = bernsteinSecondDeriv( p, s );
        const Eigen::VectorXd t_dderiv = bernsteinSecondDeriv( q, t );

        return ( Eigen::MatrixX3d( s_evals.size() * t_evals.size(), 3 ) << Eigen::kroneckerProduct( t_evals, s_dderiv ),
                 Eigen::kroneckerProduct( t_deriv, s_deriv ),
                 Eigen::kroneckerProduct( t_dderiv, s_evals ) )
            .finished();
    }

    Eigen::VectorXd
        bernsteinTP( const size_t p, const size_t q, const size_t r, const double s, const double t, const double u )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );
        const Eigen::VectorXd u_evals = bernstein( r, u );

        return Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_evals, s_evals ) );
    }

    Eigen::MatrixX3d bernsteinTPFirstDeriv(
        const size_t p, const size_t q, const size_t r, const double s, const double t, const double u )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );
        const Eigen::VectorXd u_evals = bernstein( r, u );
        const Eigen::VectorXd s_deriv = bernsteinFirstDeriv( p, s );
        const Eigen::VectorXd t_deriv = bernsteinFirstDeriv( q, t );
        const Eigen::VectorXd u_deriv = bernsteinFirstDeriv( r, u );

        return ( Eigen::MatrixX3d( s_evals.size() * t_evals.size() * u_evals.size(), 3 )
                     << Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_evals, s_deriv ) ),
                 Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_deriv, s_evals ) ),
                 Eigen::kroneckerProduct( u_deriv, Eigen::kroneckerProduct( t_evals, s_evals ) ) )
            .finished();
    }

    Eigen::MatrixXd bernsteinTPSecondDeriv(
        const size_t p, const size_t q, const size_t r, const double s, const double t, const double u )
    {
        const Eigen::VectorXd s_evals = bernstein( p, s );
        const Eigen::VectorXd t_evals = bernstein( q, t );
        const Eigen::VectorXd u_evals = bernstein( r, u );
        const Eigen::VectorXd s_deriv = bernsteinFirstDeriv( p, s );
        const Eigen::VectorXd t_deriv = bernsteinFirstDeriv( q, t );
        const Eigen::VectorXd u_deriv = bernsteinFirstDeriv( r, u );
        const Eigen::VectorXd s_dderiv = bernsteinSecondDeriv( p, s );
        const Eigen::VectorXd t_dderiv = bernsteinSecondDeriv( q, t );
        const Eigen::VectorXd u_dderiv = bernsteinSecondDeriv( r, u );

        return ( Eigen::MatrixXd( s_evals.size() * t_evals.size() * u_evals.size(), 6 )
                     << Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_evals, s_dderiv ) ),
                 Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_deriv, s_deriv ) ),
                 Eigen::kroneckerProduct( u_deriv, Eigen::kroneckerProduct( t_evals, s_deriv ) ),
                 Eigen::kroneckerProduct( u_evals, Eigen::kroneckerProduct( t_dderiv, s_evals ) ),
                 Eigen::kroneckerProduct( u_deriv, Eigen::kroneckerProduct( t_deriv, s_evals ) ),
                 Eigen::kroneckerProduct( u_dderiv, Eigen::kroneckerProduct( t_evals, s_evals ) ) )
            .finished();
    }

    enum class VectorEvalType
    {
        DivConf,
        CurlConf
    };

    Eigen::MatrixX2d bernsteinVectorConforming( const size_t p, const size_t q, const double s, const double t, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::VectorXd& vec1_evals, const Eigen::VectorXd& vec2_evals ) {
            return ( Eigen::MatrixX2d( vec1_evals.size() + vec2_evals.size(), 2 ) << vec1_evals,
                     Eigen::VectorXd::Zero( vec1_evals.size() ),
                     Eigen::VectorXd::Zero( vec2_evals.size() ),
                     vec2_evals )
                .finished();
        };

        if( type == VectorEvalType::DivConf )
        {
            const Eigen::VectorXd vec1_evals = bernsteinTP( p, q - 1, s, t );
            const Eigen::VectorXd vec2_evals = bernsteinTP( p - 1, q, s, t );

            return arrange_components( vec1_evals, vec2_evals );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::VectorXd vec1_evals = bernsteinTP( p - 1, q, s, t );
            const Eigen::VectorXd vec2_evals = bernsteinTP( p, q - 1, s, t );

            return arrange_components( vec1_evals, vec2_evals );
        }
        return Eigen::MatrixX2d::Zero( 0, 2 );
    }

    Eigen::MatrixX4d bernsteinVectorConformingFirstDeriv( const size_t p, const size_t q, const double s, const double t, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::MatrixX2d& vec1_deriv, const Eigen::MatrixX2d& vec2_deriv ) {
            return ( Eigen::MatrixX4d( vec1_deriv.rows() + vec2_deriv.rows(), 4 ) << vec1_deriv.col( 0 ),
                     Eigen::VectorXd::Zero( vec1_deriv.rows() ),
                     vec1_deriv.col( 1 ),
                     Eigen::VectorXd::Zero( vec1_deriv.rows() ),
                     Eigen::VectorXd::Zero( vec2_deriv.rows() ),
                     vec2_deriv.col( 0 ),
                     Eigen::VectorXd::Zero( vec2_deriv.rows() ),
                     vec2_deriv.col( 1 ) )
                .finished();
        };
        if( type == VectorEvalType::DivConf )
        {
            const Eigen::MatrixX2d vec1_deriv = bernsteinTPFirstDeriv( p, q - 1, s, t );
            const Eigen::MatrixX2d vec2_deriv = bernsteinTPFirstDeriv( p - 1, q, s, t );

            return arrange_components( vec1_deriv, vec2_deriv );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::MatrixX2d vec1_deriv = bernsteinTPFirstDeriv( p - 1, q, s, t );
            const Eigen::MatrixX2d vec2_deriv = bernsteinTPFirstDeriv( p, q - 1, s, t );

            return arrange_components( vec1_deriv, vec2_deriv );
        }
        return Eigen::MatrixX4d::Zero( 0, 4 );
    }

    Eigen::MatrixXd bernsteinVectorConformingSecondDeriv( const size_t p, const size_t q, const double s, const double t, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::MatrixXd& vec1_deriv, const Eigen::MatrixXd& vec2_deriv ) {
            return ( Eigen::MatrixXd( vec1_deriv.rows() + vec2_deriv.rows(), 6 ) << vec1_deriv.col( 0 ),
                    Eigen::VectorXd::Zero( vec1_deriv.rows() ),
                    vec1_deriv.col( 1 ),
                    Eigen::VectorXd::Zero( vec1_deriv.rows() ),
                    vec1_deriv.col( 2 ),
                    Eigen::VectorXd::Zero( vec1_deriv.rows() ),
                    Eigen::VectorXd::Zero( vec2_deriv.rows() ),
                    vec2_deriv.col( 0 ),
                    Eigen::VectorXd::Zero( vec2_deriv.rows() ),
                    vec2_deriv.col( 1 ),
                    Eigen::VectorXd::Zero( vec2_deriv.rows() ),
                    vec2_deriv.col( 2 ) )
                .finished();
        };

        if( type == VectorEvalType::DivConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPSecondDeriv( p, q - 1, s, t );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPSecondDeriv( p - 1, q, s, t );

            return arrange_components( vec1_deriv, vec2_deriv );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPSecondDeriv( p - 1, q, s, t );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPSecondDeriv( p, q - 1, s, t );

            return arrange_components( vec1_deriv, vec2_deriv );
        }
        return Eigen::MatrixXd::Zero( 0, 6 );
    }

    Eigen::MatrixX3d bernsteinVectorConforming( const size_t p, const size_t q, const size_t r, const double s, const double t, const double u, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::VectorXd& vec1_evals, const Eigen::VectorXd& vec2_evals, const Eigen::VectorXd& vec3_evals ) {
            return ( Eigen::MatrixX3d( vec1_evals.size() + vec2_evals.size() + vec3_evals.size(), 3 ) << vec1_evals,
                     Eigen::VectorXd::Zero( vec1_evals.size() ),
                     Eigen::VectorXd::Zero( vec1_evals.size() ),
                     Eigen::VectorXd::Zero( vec2_evals.size() ),
                     vec2_evals,
                     Eigen::VectorXd::Zero( vec2_evals.size() ),
                     Eigen::VectorXd::Zero( vec3_evals.size() ),
                     Eigen::VectorXd::Zero( vec3_evals.size() ),
                     vec3_evals )
                .finished();
        };

        if( type == VectorEvalType::DivConf )
        {
            const Eigen::VectorXd vec1_evals = bernsteinTP( p, q - 1, r - 1, s, t, u );
            const Eigen::VectorXd vec2_evals = bernsteinTP( p - 1, q, r - 1, s, t, u );
            const Eigen::VectorXd vec3_evals = bernsteinTP( p - 1, q - 1, r, s, t, u );

            return arrange_components( vec1_evals, vec2_evals, vec3_evals );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::VectorXd vec1_evals = bernsteinTP( p - 1, q, r, s, t, u );
            const Eigen::VectorXd vec2_evals = bernsteinTP( p, q - 1, r, s, t, u );
            const Eigen::VectorXd vec3_evals = bernsteinTP( p, q, r - 1, s, t, u );

            return arrange_components( vec1_evals, vec2_evals, vec3_evals );
        }
        return Eigen::MatrixX3d::Zero( 0, 3 );
    }

    Eigen::MatrixXd bernsteinVectorConformingFirstDeriv( const size_t p, const size_t q, const size_t r, const double s, const double t, const double u, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::MatrixX3d& vec1_deriv,
                                            const Eigen::MatrixX3d& vec2_deriv,
                                            const Eigen::MatrixX3d& vec3_deriv ) {
            return ( Eigen::MatrixXd( vec1_deriv.rows() + vec2_deriv.rows() + vec3_deriv.rows(), 9 )
                         << vec1_deriv.col( 0 ),
                            Eigen::MatrixXd::Zero( vec1_deriv.rows(), 2 ),
                            vec1_deriv.col( 1 ),
                            Eigen::MatrixXd::Zero( vec1_deriv.rows(), 2 ),
                            vec1_deriv.col( 2 ),
                            Eigen::MatrixXd::Zero( vec1_deriv.rows(), 2 ),
                            Eigen::MatrixXd::Zero( vec2_deriv.rows(), 1 ),
                            vec2_deriv.col( 0 ),
                            Eigen::MatrixXd::Zero( vec2_deriv.rows(), 2 ),
                            vec2_deriv.col( 1 ),
                            Eigen::MatrixXd::Zero( vec2_deriv.rows(), 2 ),
                            vec2_deriv.col( 2 ),
                            Eigen::MatrixXd::Zero( vec2_deriv.rows(), 1 ),
                            Eigen::MatrixXd::Zero( vec3_deriv.rows(), 2 ),
                            vec3_deriv.col( 0 ),
                            Eigen::MatrixXd::Zero( vec3_deriv.rows(), 2 ),
                            vec3_deriv.col( 1 ),
                            Eigen::MatrixXd::Zero( vec3_deriv.rows(), 2 ),
                            vec3_deriv.col( 2 )
                    ).finished();
        };

        if( type == VectorEvalType::DivConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPFirstDeriv( p, q - 1, r - 1, s, t, u );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPFirstDeriv( p - 1, q, r - 1, s, t, u );
            const Eigen::MatrixX3d vec3_deriv = bernsteinTPFirstDeriv( p - 1, q - 1, r, s, t, u );

            return arrange_components( vec1_deriv, vec2_deriv, vec3_deriv );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPFirstDeriv( p - 1, q, r, s, t, u );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPFirstDeriv( p, q - 1, r, s, t, u );
            const Eigen::MatrixX3d vec3_deriv = bernsteinTPFirstDeriv( p, q, r - 1, s, t, u );

            return arrange_components( vec1_deriv, vec2_deriv, vec3_deriv );
        }
        return Eigen::MatrixXd::Zero( 0, 9 );
    }

    Eigen::MatrixXd bernsteinVectorConformingSecondDeriv( const size_t p, const size_t q, const size_t r, const double s, const double t, const double u, const VectorEvalType type )
    {
        const auto arrange_components = []( const Eigen::MatrixXd& vec1_deriv,
                                            const Eigen::MatrixXd& vec2_deriv,
                                            const Eigen::MatrixXd& vec3_deriv ) {
            Eigen::MatrixXd result =
                Eigen::MatrixXd::Zero( vec1_deriv.rows() + vec2_deriv.rows() + vec3_deriv.rows(), 18 );
            result.block( 0, 0, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 0 );
            result.block( 0, 3, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 1 );
            result.block( 0, 6, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 2 );
            result.block( 0, 9, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 3 );
            result.block( 0, 12, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 4 );
            result.block( 0, 15, vec1_deriv.rows(), 1 ) = vec1_deriv.col( 5 );
            result.block( vec1_deriv.rows(), 1, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 0 );
            result.block( vec1_deriv.rows(), 4, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 1 );
            result.block( vec1_deriv.rows(), 7, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 2 );
            result.block( vec1_deriv.rows(), 10, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 3 );
            result.block( vec1_deriv.rows(), 13, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 4 );
            result.block( vec1_deriv.rows(), 16, vec2_deriv.rows(), 1 ) = vec2_deriv.col( 5 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 2, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 0 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 5, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 1 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 8, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 2 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 11, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 3 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 14, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 4 );
            result.block( vec1_deriv.rows() + vec2_deriv.rows(), 17, vec3_deriv.rows(), 1 ) = vec3_deriv.col( 5 );
            return result;
        };

        if( type == VectorEvalType::DivConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPSecondDeriv( p, q - 1, r - 1, s, t, u );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPSecondDeriv( p - 1, q, r - 1, s, t, u );
            const Eigen::MatrixX3d vec3_deriv = bernsteinTPSecondDeriv( p - 1, q - 1, r, s, t, u );

            return arrange_components( vec1_deriv, vec2_deriv, vec3_deriv );
        }
        else if( type == VectorEvalType::CurlConf )
        {
            const Eigen::MatrixX3d vec1_deriv = bernsteinTPSecondDeriv( p - 1, q, r, s, t, u );
            const Eigen::MatrixX3d vec2_deriv = bernsteinTPSecondDeriv( p, q - 1, r, s, t, u );
            const Eigen::MatrixX3d vec3_deriv = bernsteinTPSecondDeriv( p, q, r - 1, s, t, u );

            return arrange_components( vec1_deriv, vec2_deriv, vec3_deriv );
        }
        return Eigen::MatrixXd::Zero( 0, 18 );
    }

    inline Eigen::Index numCols( const size_t dim, const size_t n_deriv )
    {
        switch( dim )
        {
            case 1: return n_deriv + 1;
            case 2: return ( n_deriv + 1 ) * ( n_deriv + 2 ) / 2;
            case 3: return ( n_deriv + 1 ) * ( n_deriv + 2 ) * ( n_deriv + 3 ) / 6;
            default: throw std::runtime_error( "Unsupported dimension: " + std::to_string( dim ) );
        }
    }

    ParentBasisEval::ParentBasisEval( const basis::ParentBasis& pb,
                                      const param::ParentPoint& point,
                                      const size_t n_derivatives )
    {
        if( not isCartesian( pb.mParentDomain ) ) throw std::runtime_error( "Evals only implemented for TP cells" );

        const size_t param_dim = dim( pb.mParentDomain );
        const size_t n_components = numVectorComponents( pb );
        const SmallVector<size_t, 3> degs = degrees( pb );
        mEvals = Eigen::MatrixXd( numFunctions( pb ), numCols( param_dim, n_derivatives ) * n_components );

        if( n_components == 1 )
        {
            switch( param_dim )
            {
                case 1:
                    mEvals.col( 0 ) = bernstein( degs.at( 0 ), point.mPoint( 0 ) );
                    if( n_derivatives > 0 )
                    {
                        mEvals.col( 1 ) = bernsteinFirstDeriv( degs.at( 0 ), point.mPoint( 0 ) );
                        if( n_derivatives > 1 ) mEvals.col( 2 ) = bernsteinSecondDeriv( degs.at( 0 ), point.mPoint( 0 ) );
                    }
                    break;
                case 2:
                    mEvals.col( 0 ) = bernsteinTP( degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ) );
                    if( n_derivatives > 0 )
                    {
                        mEvals.middleCols<2>( 1 ) =
                            bernsteinTPFirstDeriv( degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ) );
                        if( n_derivatives > 1 )
                            mEvals.middleCols<3>( 3 ) =
                                bernsteinTPSecondDeriv( degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ) );
                    }
                    break;
                case 3:
                    mEvals.col( 0 ) = bernsteinTP(
                        degs.at( 0 ), degs.at( 1 ), degs.at( 2 ), point.mPoint( 0 ), point.mPoint( 1 ), point.mPoint( 2 ) );
                    if( n_derivatives > 0 )
                    {
                        mEvals.middleCols<3>( 1 ) = bernsteinTPFirstDeriv( degs.at( 0 ),
                                                                        degs.at( 1 ),
                                                                        degs.at( 2 ),
                                                                        point.mPoint( 0 ),
                                                                        point.mPoint( 1 ),
                                                                        point.mPoint( 2 ) );
                        if( n_derivatives > 1 )
                            mEvals.middleCols<6>( 4 ) = bernsteinTPSecondDeriv( degs.at( 0 ),
                                                                                degs.at( 1 ),
                                                                                degs.at( 2 ),
                                                                                point.mPoint( 0 ),
                                                                                point.mPoint( 1 ),
                                                                                point.mPoint( 2 ) );
                    }
                    break;
                default: throw std::runtime_error( "Unsupported dimension: " + std::to_string( param_dim ) );
            }
        }
        else
        {
            if( pb.mBasisGroups.size() != 1 )
                throw std::runtime_error( "Vector-valued splines cannot be tensor producted." );

            const VectorEvalType eval_type =
                pb.mBasisGroups.at( 0 ).type == basis::BasisType::DivConformingBernstein
                    ? VectorEvalType::DivConf
                    : VectorEvalType::CurlConf;
            switch( param_dim )
            {
                case 2:
                    mEvals.leftCols( 2 ) = bernsteinVectorConforming(
                        degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ), eval_type );
                    if( n_derivatives > 0 )
                    {
                        mEvals.middleCols<4>( 2 ) = bernsteinVectorConformingFirstDeriv(
                            degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ), eval_type );
                        if( n_derivatives > 1 )
                            mEvals.middleCols<6>( 6 ) = bernsteinVectorConformingSecondDeriv(
                                degs.at( 0 ), degs.at( 1 ), point.mPoint( 0 ), point.mPoint( 1 ), eval_type );
                    }
                    break;
                case 3:
                    mEvals.leftCols( 3 ) = bernsteinVectorConforming( degs.at( 0 ),
                                                                      degs.at( 1 ),
                                                                      degs.at( 2 ),
                                                                      point.mPoint( 0 ),
                                                                      point.mPoint( 1 ),
                                                                      point.mPoint( 2 ),
                                                                      eval_type );
                    if( n_derivatives > 0 )
                    {
                        mEvals.middleCols<9>( 3 ) = bernsteinVectorConformingFirstDeriv( degs.at( 0 ),
                                                                                         degs.at( 1 ),
                                                                                         degs.at( 2 ),
                                                                                         point.mPoint( 0 ),
                                                                                         point.mPoint( 1 ),
                                                                                         point.mPoint( 2 ),
                                                                                         eval_type );
                        if( n_derivatives > 1 )
                            mEvals.middleCols<18>( 12 ) = bernsteinVectorConformingSecondDeriv( degs.at( 0 ),
                                                                                                degs.at( 1 ),
                                                                                                degs.at( 2 ),
                                                                                                point.mPoint( 0 ),
                                                                                                point.mPoint( 1 ),
                                                                                                point.mPoint( 2 ),
                                                                                                eval_type );
                    }
                    break;
                default: throw std::runtime_error( "Unsupported dimension for vector-valued splines: " + std::to_string( param_dim ) );
            }


        }
    }
} // namespace eval