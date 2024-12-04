#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <numbers>
#include <iostream>
#include <Eigen/LU>

namespace py = pybind11;
using namespace py::literals;
using std::numbers::pi;

Eigen::MatrixXd eigenfunction_series(
    const Eigen::MatrixXd& x,
    const Eigen::MatrixXd& y,
    const size_t n_terms,
    const std::function<double( const size_t, const size_t )>&
        func,
    const double multiplier )
{
    Eigen::MatrixXd u( x.rows(), x.cols() );
    
    for( size_t n = 1; n <= n_terms; n++ )
        for( size_t m = 1; m <= n_terms; m++ )
            u += func( m, n ) * ( Eigen::sin( ( m * pi * x.array() ) / 9 ) * Eigen::sin( ( n * pi * y.array() ) / 3 ) ).matrix();

    u *= multiplier;
    return u;
}

Eigen::MatrixXd eigenfunction_series(
    const Eigen::MatrixXd& x,
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& z,
    const size_t n_terms,
    const std::function<Eigen::MatrixXd( const Eigen::MatrixXd&,
    const Eigen::MatrixXd&,
    const Eigen::MatrixXd&, const size_t, const size_t )>& func )
{
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero( x.rows(), x.cols() );
    // std::cout << "0 0 " << u.transpose() << std::endl;
    
    for( size_t n = 1; n <= n_terms; n++ )
        for( size_t m = 0; m <= n_terms; m++ )
        {
            u += func( x, y, z, m, n );
            if( Eigen::isnan( u.array() ).any() )
            {
                std::cout << n << " " << m << " has a nan\n";
                return u;
            }
            if( Eigen::isinf( u.array() ).any() )
            {
                std::cout << n << " " << m << " has a nan\n";
                return u;
            }
            // std::cout << m << " " << n << " " << u.transpose() << std::endl;
        }

    return u;
}

double u11( const size_t m, const size_t n )
{
    return ( 1.0 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( -162 + 49 * pow( m, 2 ) * pow( pi, 2 ) ) * cos( ( m * pi ) / 9 ) +
             18 * ( 9 * cos( 8 * m * pi / 9 ) + 7 * m * pi * sin( ( m * pi ) / 9 ) ) ) *
           sin( n * pi / 3 );
}

double u12( const size_t m, const size_t n )
{

    return -( 36 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( 7 * m * pi * cos( 7 * m * pi / 18 ) - 18 * sin( 7 * m * pi / 18 ) ) *
           sin( m * pi / 2 ) * sin( n * pi / 3 );
}

double u13( const size_t m, const size_t n )
{
    return -( 1.0 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( 162 * cos( m * pi / 9 ) +
             ( -162 + 49 * pow( m, 2 ) * pow( pi, 2 ) ) * cos( 8 * m * pi / 9 ) -
             126 * m * pi * sin( 8 * m * pi / 9 ) ) *
           sin( n * pi / 3 );
}

double u21( const size_t m, const size_t n )
{
    return ( 1.0 / ( 3 * pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( pow( n, 2 ) * pow( pi, 2 ) - 18 ) * cos( n * pi / 3 ) +
             6 * ( 3 * cos( 2 * n * pi / 3 ) + n * pi * sin( n * pi / 3 ) ) ) *
           sin( m * pi / 9 );
}

double u22( const size_t m, const size_t n )
{
    return ( -4.0 / ( pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( n * pi * cos( n * pi / 6 ) - 6 * sin( n * pi / 6 ) ) * sin( n * pi / 2 ) *
           sin( m * pi / 9 );
}

double u23( const size_t m, const size_t n )
{
    return ( -1.0 / ( 3 * pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( pow( n, 2 ) * pow( pi, 2 ) - 18 ) * cos( 2 * n * pi / 3 ) +
             6 * ( 3 * cos( n * pi / 3 ) - n * pi * sin( 2 * n * pi / 3 ) ) ) *
           sin( m * pi / 9 );
}

double u31( const size_t m, const size_t n )
{
    return ( 1.0 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( -162 + 49 * pow( m, 2 ) * pow( pi, 2 ) ) * cos( ( m * pi ) / 9 ) +
             18 * ( 9 * cos( 8 * m * pi / 9 ) + 7 * m * pi * sin( ( m * pi ) / 9 ) ) ) *
           sin( 2 * n * pi / 3 );
}

double u32( const size_t m, const size_t n )
{
    return -( 36 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( 7 * m * pi * cos( 7 * m * pi / 18 ) - 18 * sin( 7 * m * pi / 18 ) ) *
           sin( m * pi / 2 ) * sin( 2 * n * pi / 3 );
}

double u33( const size_t m, const size_t n )
{
    return -( 1.0 / ( 49 * pow( m, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( 162 * cos( m * pi / 9 ) +
             ( -162 + 49 * pow( m, 2 ) * pow( pi, 2 ) ) * cos( 8 * m * pi / 9 ) -
             126 * m * pi * sin( 8 * m * pi / 9 ) ) *
           sin( 2 * n * pi / 3 );
}

double u41( const size_t m, const size_t n )
{
    return ( 1.0 / ( 3 * pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( pow( n, 2 ) * pow( pi, 2 ) - 18 ) * cos( n * pi / 3 ) +
             6 * ( 3 * cos( 2 * n * pi / 3 ) + n * pi * sin( n * pi / 3 ) ) ) *
           sin( 8 * m * pi / 9 );
}

double u42( const size_t m, const size_t n )
{
    return ( -4.0 / ( pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( n * pi * cos( n * pi / 6 ) - 6 * sin( n * pi / 6 ) ) * sin( n * pi / 2 ) *
           sin( 8 * m * pi / 9 );
}

double u43( const size_t m, const size_t n )
{
    return ( -1.0 / ( 3 * pow( n, 3 ) * ( pow( m, 2 ) + 9 * pow( n, 2 ) ) * pow( pi, 3 ) ) ) *
           ( ( pow( n, 2 ) * pow( pi, 2 ) - 18 ) * cos( 2 * n * pi / 3 ) +
             6 * ( 3 * cos( n * pi / 3 ) - n * pi * sin( 2 * n * pi / 3 ) ) ) *
           sin( 8 * m * pi / 9 );
}

std::function<double( const size_t, const size_t )>
    G( const double xi, const double eta )
{
    return
        [xi, eta]( const size_t m, const size_t n ) -> double {
            return sin( m * pi * xi / 9 ) * sin( n * pi * eta / 3 ) / ( 9 * pow( m, 2 ) + 81 * pow( n, 2 ) );
        };
}

Eigen::MatrixXd gamma( const Eigen::MatrixXd& vals, const bool at_one, const bool is_x, const size_t m, const size_t n )
{
    const double xval = at_one ? 1.0 : 0.0;
    const int factor = is_x ? 1.0 : -1.0;
    const size_t factor2 = not is_x and m == 0 ? 2 : 1;
    const double beta = sqrt( pow( m, 2 ) + pow( is_x ? ( 2 * n - 1 ) : n, 2 ) );
    return 4 * ( Eigen::exp( beta * pi * ( vals.array() - 2.0 + xval ) ) + factor * Eigen::exp( -beta * pi * ( vals.array() + xval ) ) ).matrix() *
            ( 1.0 / tanh( beta * pi ) + 1 ) /
            ( beta * beta * pi * pi * factor2 );
}

Eigen::MatrixXd z( const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::MatrixXd& Z, const size_t m, const size_t n, const double zz )
{
    return gamma( Z, true, false, m, n ).array() * Eigen::cos( m * pi * X.array() ) * Eigen::sin( n * pi * Y.array() ) * zz;
}

double z1( const size_t m, const size_t n )
{
    if( m == 0 ) return ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / 3;
    else return 2 * ( m * pi - sin( m * pi ) ) * ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( m * pi, 3 );
}

double z2( const size_t m, const size_t n )
{
    if( m == 0 ) return ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / 3;
    else return -2 * ( m * pi + m * pi * cos( m * pi ) - 2 * sin( m * pi ) ) * ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( m * pi, 3 );
}

double z3( const size_t m, const size_t n )
{
    if( m == 0 ) return ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / 3;
    else return ( 2 * m * pi * cos( m * pi ) + ( m * m * pi * pi - 2 ) * sin( m * pi ) ) * ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( m * pi, 3 );
}

Eigen::MatrixXd x0( const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::MatrixXd& Z, const size_t n, const size_t k, const double xx )
{
    return gamma( X, false, true, n, k ).array() * Eigen::sin( n * pi * Y.array() ) * Eigen::sin( ( 2 * k - 1 ) * pi * Z.array() / 2 ) * xx;
}

Eigen::MatrixXd x1( const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::MatrixXd& Z, const size_t n, const size_t k, const double xx )
{
    return gamma( X, true, true, n, k ).array() * Eigen::sin( n * pi * Y.array() ) * Eigen::sin( ( 2 * k - 1 ) * pi * Z.array() / 2 ) * xx;
}

double xx1( size_t n, size_t k )
{
    if( n == 0 ) return 0;
    const size_t kk = 2 * k - 1;
    return 2 * ( 6 * kk * pi * cos( ( k - 2 ) * pi / 3 ) + ( kk * kk * pi * pi - 18 ) * cos( kk * pi / 6 ) + 18 * sin( k * pi ) ) *
        ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( pi * kk, 3 );
}

double xx2( size_t n, size_t k )
{
    if( n == 0 ) return 0;
    const size_t kk = 2 * k - 1;
    return 24 * sin( kk * pi / 3 ) * ( 6 * cos( ( k - 2 ) * pi / 3 ) - kk * pi * cos( kk * pi / 6 ) ) *
        ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( pi * kk, 3 );
}

double xx3( size_t n, size_t k )
{
    if( n == 0 ) return 0;
    const size_t kk = 2 * k - 1;
    return -2 * ( 6 * kk * pi * cos( k * pi ) + 18 * cos( kk * pi / 6 ) + ( kk * kk * pi * pi - 18 ) * sin( k * pi ) ) *
        ( sin( n * pi / 3 ) + sin( 2 * n * pi / 3 ) ) / pow( pi * kk, 3 );
}

double xy1( size_t n, size_t k )
{
    if( n == 0 ) return 0;
    return 2 * cos( ( k - 2 ) * pi / 3 ) * sin( n * pi / 2 ) * ( 6 * n * pi * cos( n * pi / 6 ) + ( n * n * pi * pi - 36 ) * sin( n * pi / 6 ) ) / pow( n * pi, 3 );
}

double xy2( size_t n, size_t k )
{
    if( n == 0 ) return 0;
    return -12 * cos( ( k - 2 ) * pi / 3 ) * sin( n * pi / 2 ) * ( n * pi * cos( n * pi / 6 ) - 6 * sin( n * pi / 6 ) ) / pow( n * pi, 3 );
}

double xmid( size_t n, size_t k )
{
    return sin( ( 2 * k - 1 ) * pi / 4 ) * sin( n * pi / 2 );
}

constexpr double Lx = 9;
constexpr double Ly = 3;
const double mult = pow( Lx, 3 ) * pow( Ly, 3 ) / 4 / pow( pi, 2 );

PYBIND11_MODULE( eigenfunction, m )
{
    m.doc() = "Plugin for performing faster calculations of eigenfunction expansions";

    m.def( "leastSquaresFit",
           []( const size_t n_terms, const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& values ) {
                Eigen::MatrixXd A( X.rows(), 14 );
                
                A.col( 0 ) = eigenfunction_series(X, Y, n_terms, u11, mult);
                A.col( 1 ) = eigenfunction_series(X, Y, n_terms, u12, mult);
                A.col( 2 ) = eigenfunction_series(X, Y, n_terms, u13, mult);
                A.col( 3 ) = eigenfunction_series(X, Y, n_terms, u21, mult);
                A.col( 4 ) = eigenfunction_series(X, Y, n_terms, u22, mult);
                A.col( 5 ) = eigenfunction_series(X, Y, n_terms, u23, mult);
                A.col( 6 ) = eigenfunction_series(X, Y, n_terms, u31, mult);
                A.col( 7 ) = eigenfunction_series(X, Y, n_terms, u32, mult);
                A.col( 8 ) = eigenfunction_series(X, Y, n_terms, u33, mult);
                A.col( 9 ) = eigenfunction_series(X, Y, n_terms, u41, mult);
                A.col( 10 ) = eigenfunction_series(X, Y, n_terms, u42, mult);
                A.col( 11 ) = eigenfunction_series(X, Y, n_terms, u43, mult);
                A.col( 12 ) = eigenfunction_series(X, Y, n_terms, G( 1.5, 1.5 ), mult);
                A.col( 13 ) = eigenfunction_series(X, Y, n_terms, G( 7.5, 1.5 ), mult);
                
                std::cout << "Running solve:\n";
                std::cout << A << std::endl << std::endl;
                std::cout << values.transpose() << std::endl;

                Eigen::PartialPivLU<Eigen::MatrixXd> sol( A );

                const Eigen::VectorXd c = sol.solve( values );
                return c;
           },
           "Computes the level sets of the sweep parameterization and traces the parameterization, outputting the "
           "results to .vtu files.",
           "n_terms"_a,
           "X"_a,
           "Y"_a,
           "values"_a );

    m.def( 
        "leastSquaresFitVolume",
        []( const size_t n_terms, const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& Z, const Eigen::VectorXd& values ){
            Eigen::MatrixXd A( X.rows(), 7 );

            A.col( 0 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return z( X, Y, Z, m, n, z1( m, n ) + z3( m, n ) ) + x0( X, Y, Z, m, n, xx3( m, n ) ) + x1( X, Y, Z, m, n, xx3( m, n ) );
            } );

            A.col( 1 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return z( X, Y, Z, m, n, z2( m, n ) );
            } );

            A.col( 2 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return x0( X, Y, Z, m, n, xx2( m, n ) ) + x1( X, Y, Z, m, n, xx2( m, n ) );
            } );

            A.col( 3 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return x0( X, Y, Z, m, n, xx1( m, n ) ) + x1( X, Y, Z, m, n, xx1( m, n ) );
            } );

            A.col( 4 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return x0( X, Y, Z, m, n, xy1( m, n ) ) + x1( X, Y, Z, m, n, xy1( m, n ) );
            } );

            A.col( 5 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return x0( X, Y, Z, m, n, xy2( m, n ) ) + x1( X, Y, Z, m, n, xy2( m, n ) );
            } );

            A.col( 6 ) = eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                return x0( X, Y, Z, m, n, xmid( m, n ) ) + x1( X, Y, Z, m, n, xmid( m, n ) );
            } );

            std::cout << "Running solve:\n";
            std::cout << A << std::endl << std::endl;
            std::cout << values.transpose() << std::endl;

            Eigen::FullPivLU<Eigen::MatrixXd> sol( A );

            const Eigen::VectorXd c = sol.solve( values );

            std::cout << ( A*c ).transpose() << std::endl;
            return c;
            return Eigen::VectorXd();
        } );

    m.def(
        "evaluateFit",
        []( const size_t n_terms, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& c ) {
            const auto solution_terms =
                [&c]( const size_t m, const size_t n ) {
                    return c[0] * u11( m, n ) + c( 1 ) * u12( m, n ) + c( 2 ) * u13( m, n ) +
                           c( 3 ) * u21( m, n ) + c( 4 ) * u22( m, n ) + c( 5 ) * u23( m, n ) +
                           c( 6 ) * u31( m, n ) + c( 7 ) * u32( m, n ) + c( 8 ) * u33( m, n ) +
                           c( 9 ) * u41( m, n ) + c( 10 ) * u42( m, n ) + c( 11 ) * u43( m, n ) +
                           c( 12 ) * G( 1.5, 1.5 )( m, n ) + c( 13 ) * G( 7.5, 1.5 )( m, n );
                };
            
            return eigenfunction_series( X, Y, n_terms, solution_terms, mult );
        },
        "Computes the level sets of the sweep parameterization and traces the parameterization, outputting the "
        "results to .vtu files.",
        "n_terms"_a,
        "X"_a,
        "Y"_a,
        "c"_a );

    m.def( 
        "evaluateFitVolume",
        []( const size_t n_terms, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::MatrixXd& Z, const Eigen::VectorXd& c ){
            return eigenfunction_series( X, Y, Z, n_terms, [&]( const auto& X, const auto& Y, const auto& Z, const size_t m, const size_t n ) -> Eigen::MatrixXd {
                Eigen::MatrixXd out = c( 0 ) * ( z( X, Y, Z, m, n, z1( m, n ) + z3( m, n ) ) + x0( X, Y, Z, m, n, xx3( m, n ) ) + x1( X, Y, Z, m, n, xx3( m, n ) ) );

                out += c( 1 ) * z( X, Y, Z, m, n, z2( m, n ) );
                out += c( 2 ) * ( x0( X, Y, Z, m, n, xx2( m, n ) ) + x1( X, Y, Z, m, n, xx2( m, n ) ) );
                out += c( 3 ) * ( x0( X, Y, Z, m, n, xx1( m, n ) ) + x1( X, Y, Z, m, n, xx1( m, n ) ) );
                out += c( 4 ) * ( x0( X, Y, Z, m, n, xy1( m, n ) ) + x1( X, Y, Z, m, n, xy1( m, n ) ) );
                out += c( 5 ) * ( x0( X, Y, Z, m, n, xy2( m, n ) ) + x1( X, Y, Z, m, n, xy2( m, n ) ) );
                out += c( 6 ) * ( x0( X, Y, Z, m, n, xmid( m, n ) ) + x1( X, Y, Z, m, n, xmid( m, n ) ) );

                return out;
            } );
        } );
}