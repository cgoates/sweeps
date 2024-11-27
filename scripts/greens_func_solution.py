import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from numpy import pi
import sys
from pathlib import Path
path_to_api = Path(__file__).parent.parent.parent / "build" / "src" / "api"
sys.path.insert(0, "/Users/caleb/sweeps/build/src/api")
import eigenfunction

np.set_printoptions(linewidth=270)

def eigenfunction_series(x, y, n_terms, func, multiplier):
    u = np.zeros_like(x, dtype=np.float64)
    
    for n in range(1, n_terms+1):
        for m in range(1, n_terms+1):
            # Eigenfunction φ_nm(x, y)
            u += func( x, y, m, n )
    
    u *= multiplier
    return u

def constant_along_interior_bdry( x, y, m, n ):
    Lx = 9
    Ly = 3
    return np.sin( m * np.pi * x / Lx ) * np.sin( n * np.pi * y / Ly ) / ( m**2 * Ly**2 + n**2 * Lx**2 )\
                * ( ( np.sin(n*np.pi/3) + np.sin(2*n*np.pi/3) ) * 9 / m / np.pi * ( np.cos( m * np.pi/9 ) - np.cos( 8*m*np.pi/9 ) )\
                    + ( np.sin(m*np.pi/9) + np.sin(8*m*np.pi/9) ) * 9 / m / np.pi * ( np.cos( n * np.pi/3 ) - np.cos( 2*n*np.pi/3 ) ) )

def u11( x, y, m, n ):
    return ( 1.0 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( -162 + 49 * m**2 * pi**2 ) * np.cos( ( m * pi ) / 9 ) + 18 * ( 9 * np.cos(8*m*pi/9) + 7 * m * pi * np.sin( ( m * pi ) / 9 ) ) ) * \
            np.sin( n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u12( x, y, m, n ):
    return -( 36 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( 7 * m * pi * np.cos( 7 * m * pi / 18 ) - 18 * np.sin( 7 * m * pi / 18 ) ) * np.sin( m * pi / 2 ) * \
        np.sin( n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u13( x, y, m, n ):
    return -( 1.0 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( 162 * np.cos( m * pi / 9 ) + ( -162 + 49 * m**2 * pi**2 ) * np.cos( 8 * m * pi / 9 ) - 126 * m * pi * np.sin( 8 * m * pi / 9 ) ) * \
        np.sin( n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u21( x, y, m, n ):
    return ( 1.0 / ( 3 * n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( n**2 * pi**2 - 18 ) * np.cos( n * pi / 3 ) + 6 * ( 3 * np.cos( 2 * n * pi / 3 ) + n * pi * np.sin( n * pi / 3 ) ) ) * \
        np.sin( m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )

def u22( x, y, m, n ):
    return ( - 4.0 / ( n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( n * pi * np.cos( n * pi / 6 ) - 6 * np.sin( n * pi / 6 ) ) * np.sin( n * pi / 2 ) * \
        np.sin( m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )

def u23( x, y, m, n ):
    return ( - 1.0 / ( 3 * n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( n**2 * pi**2 - 18 ) * np.cos( 2 * n * pi / 3 ) + 6 * ( 3 * np.cos( n * pi / 3 ) - n * pi * np.sin( 2 * n * pi / 3 ) ) ) * \
        np.sin( m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )

def u31( x, y, m, n ):
    return ( 1.0 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( -162 + 49 * m**2 * pi**2 ) * np.cos( ( m * pi ) / 9 ) + 18 * ( 9 * np.cos(8*m*pi/9) + 7 * m * pi * np.sin( ( m * pi ) / 9 ) ) ) * \
            np.sin( 2 * n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u32( x, y, m, n ):
    return -( 36 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( 7 * m * pi * np.cos( 7 * m * pi / 18 ) - 18 * np.sin( 7 * m * pi / 18 ) ) * np.sin( m * pi / 2 ) * \
        np.sin( 2 * n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u33( x, y, m, n ):
    return -( 1.0 / ( 49 * m**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( 162 * np.cos( m * pi / 9 ) + ( -162 + 49 * m**2 * pi**2 ) * np.cos( 8 * m * pi / 9 ) - 126 * m * pi * np.sin( 8 * m * pi / 9 ) ) * \
        np.sin( 2 * n * pi / 3 ) * np.sin( ( m * pi * x ) / 9 ) * np.sin( ( n * pi * y ) / 3 )

def u41( x, y, m, n ):
    return ( 1.0 / ( 3 * n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( n**2 * pi**2 - 18 ) * np.cos( n * pi / 3 ) + 6 * ( 3 * np.cos( 2 * n * pi / 3 ) + n * pi * np.sin( n * pi / 3 ) ) ) * \
        np.sin( 8 * m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )

def u42( x, y, m, n ):
    return ( - 4.0 / ( n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( n * pi * np.cos( n * pi / 6 ) - 6 * np.sin( n * pi / 6 ) ) * np.sin( n * pi / 2 ) * \
        np.sin( 8 * m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )

def u43( x, y, m, n ):
    return ( - 1.0 / ( 3 * n**3 * ( m**2 + 9 * n**2 ) * pi**3 ) ) * \
        ( ( n**2 * pi**2 - 18 ) * np.cos( 2 * n * pi / 3 ) + 6 * ( 3 * np.cos( n * pi / 3 ) - n * pi * np.sin( 2 * n * pi / 3 ) ) ) * \
        np.sin( 8 * m * pi / 9 ) * np.sin( m * pi * x / 9 ) * np.sin( n * pi * y / 3 )
        
def G( xi, eta ):
    def func( x, y, m, n ):
        return np.sin( m * pi * x / 9 ) * np.sin( m * pi * xi / 9 ) * np.sin( n * pi * y / 3 ) * np.sin( n * pi * eta / 3 ) / ( 9 * m**2 + 81 * n**2 )
    
    return func

    
def evaluate_fit_solution( c, Lx, Ly, terms ):
    def fit_solution( x, y, m, n ):
        return c[0] * u11( x, y, m, n ) + \
            c[1] * u12( x, y, m, n ) + \
            c[2] * u13( x, y, m, n ) + \
            c[3] * u21( x, y, m, n ) + \
            c[4] * u22( x, y, m, n ) + \
            c[5] * u23( x, y, m, n ) + \
            c[6] * u31( x, y, m, n ) + \
            c[7] * u32( x, y, m, n ) + \
            c[8] * u33( x, y, m, n ) + \
            c[9] * u41( x, y, m, n ) + \
            c[10] * u42( x, y, m, n ) + \
            c[11] * u43( x, y, m, n ) + \
            c[12] * G( 1.5, 1.5 )( x, y, m, n ) + \
            c[13] * G( 7.5, 1.5 )( x, y, m, n )
            
    evaluate_fit_solution_along_constant_lines( fit_solution, Lx, Ly, terms )
    
def evaluate_fit_solution_along_constant_lines( fit_solution, Lx, Ly, terms ):
    X2, Y2 = np.meshgrid(np.array([1, 2, 3, 4, 4.5, 5, 6, 7, 8]), np.array([1, 1.25, 1.5, 1.75, 2]))
    u2 = eigenfunction_series(X2, Y2, terms, fit_solution, Lx**3*Ly**3/4/pi**2)
    print("Solution along line: ", u2)
    
def fit_to_constant_along_lines():
    Lx = 9
    Ly = 3
    n_terms = 1000
    # mult = Lx**3*Ly**3/4/pi**2
    X = np.array( [1, 2.5, 4.5, 6.5, 8, 1, 8, 1, 2.5, 4.5, 6.5, 8, 1.5, 7.5] )
    Y = np.array( [1, 1, 1, 1, 1, 1.5, 1.5, 2, 2, 2, 2, 2, 1.5, 1.5] )
    
    # A = np.transpose( np.vstack((
    #     eigenfunction_series(X, Y, n_terms, u11, mult),
    #     eigenfunction_series(X, Y, n_terms, u12, mult),
    #     eigenfunction_series(X, Y, n_terms, u13, mult),
    #     eigenfunction_series(X, Y, n_terms, u21, mult),
    #     eigenfunction_series(X, Y, n_terms, u22, mult),
    #     eigenfunction_series(X, Y, n_terms, u23, mult),
    #     eigenfunction_series(X, Y, n_terms, u31, mult),
    #     eigenfunction_series(X, Y, n_terms, u32, mult),
    #     eigenfunction_series(X, Y, n_terms, u33, mult),
    #     eigenfunction_series(X, Y, n_terms, u41, mult),
    #     eigenfunction_series(X, Y, n_terms, u42, mult),
    #     eigenfunction_series(X, Y, n_terms, u43, mult),
    #     eigenfunction_series(X, Y, n_terms, G( 1.5, 1.5 ), mult),
    #     eigenfunction_series(X, Y, n_terms, G( 7.5, 1.5 ), mult)
    # )) )
    
    # print(A)
    b = np.array([1.03]*12 + [4]*2)
    # c = np.linalg.lstsq(A, b )[0]
    # print( "Coefficients: ", c)
    
    c2 = eigenfunction.leastSquaresFit( n_terms, X, Y, b )
    print( "C++ Coefficients: ", c2 )
    
    X2, Y2 = np.meshgrid(np.array([1, 2, 3, 4, 4.5, 5, 6, 7, 8]), np.array([1, 1.25, 1.5, 1.75, 2]))
    u2 = eigenfunction.evaluateFit( n_terms, X2, Y2, c2 )
    print( u2 )
    
    
    # Discretize the domain (grid)
    Nx = 61  # number of grid points in x
    Ny = 21  # number of grid points in y
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Evaluate the solution u(x, y, t)
    u = eigenfunction.evaluateFit( n_terms, X, Y, c2 )

    # Plot the solution
    fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111,projection='3d')
    ax = fig.add_subplot(111)
    # cp = ax.plot_surface(X, Y, u, cmap='viridis')
    ax.axis("equal")
    divnorm = colors.TwoSlopeNorm( vmin=0, vcenter=1, vmax=4)
    cp = plt.contourf(X, Y, u, 100, cmap='RdBu_r', norm=divnorm)
    plt.contour(X, Y, u, [1.0])
    plt.colorbar(cp)
    plt.title('Eigenfunction Series Solution $u(x, y)$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


    # evaluate_fit_solution( c, Lx, Ly, n_terms )
    
def fit_to_constant_along_lines_volume():
    n_terms = 100
    X = np.array( [0]*6 + [0.5] )
    Y = np.array( [0.5]*2 + [2.0/3.0]*5 )
    Z = np.array( [ 0.5, 1.0/3, 1.0/3, 0.5, 5.0/6, 1.0, 1.0] )

    b = np.array([4.0] + [1.0]*6)
    
    c2 = eigenfunction.leastSquaresFitVolume( n_terms, X, Y, Z, b )
    print( "C++ Coefficients: ", c2 )

    # evaluate_fit_solution( c, Lx, Ly, n_terms )

fit_to_constant_along_lines_volume()