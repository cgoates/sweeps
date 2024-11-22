import numpy as np
import matplotlib.pyplot as plt
from numpy import pi

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


def evaluate_for_constant_lines():
    # Domain size
    Lx = 9
    Ly = 3

    # Discretize the domain (grid)
    Nx = 271  # number of grid points in x
    Ny = 91  # number of grid points in y
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Number of terms to consider in the series
    terms = 400

    # Evaluate the solution u(x, y, t)
    u = eigenfunction_series(X, Y, terms, constant_along_interior_bdry, Lx**3*Ly**3/np.pi**2/200)

    # Plot the solution
    fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111,projection='3d')
    # cp = ax.plot_surface(X, Y, u, cmap='viridis')
    # ax.axis("equal")
    cp = plt.contourf(X, Y, u, 100, cmap='viridis')
    plt.colorbar(cp)
    plt.title('Eigenfunction Series Solution $u(x, y)$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    
def evaluate_along_constant_lines():
    Lx = 9
    Ly = 3
    terms = 400
    X2, Y2 = np.meshgrid(np.array([1, 2, 3, 4, 5, 6, 7, 8]), np.array([1, 1.5, 2]))
    u2 = eigenfunction_series(X2, Y2, terms, constant_along_interior_bdry, Lx**3*Ly**3/np.pi**2/200)
    print(u2)
    
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

    # Discretize the domain (grid)
    Nx = 181  # number of grid points in x
    Ny = 61  # number of grid points in y
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Evaluate the solution u(x, y, t)
    u = eigenfunction_series(X, Y, terms, fit_solution, Lx**3*Ly**3/4/pi**2)

    # Plot the solution
    fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111,projection='3d')
    ax = fig.add_subplot(111)
    # cp = ax.plot_surface(X, Y, u, cmap='viridis')
    ax.axis("equal")
    cp = plt.contourf(X, Y, u, 100, cmap='RdBu_r')
    plt.colorbar(cp)
    plt.title('Eigenfunction Series Solution $u(x, y)$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    
def evaluate_fit_solution_along_constant_lines( fit_solution, Lx, Ly, terms ):
    X2, Y2 = np.meshgrid(np.array([1, 2, 3, 4, 5, 6, 7, 8]), np.array([1, 1.25, 1.5, 1.75, 2]))
    u2 = eigenfunction_series(X2, Y2, terms, fit_solution, Lx**3*Ly**3/4/pi**2)
    print("Solution along line: ", u2)
    
def fit_to_constant_along_lines():
    Lx = 9
    Ly = 3
    n_terms = 300
    mult = Lx**3*Ly**3/4/pi**2
    X = np.array( [1, 2.5, 4.5, 6.5, 8, 1, 8, 1, 2.5, 4.5, 6.5, 8, 1.5, 7.5] )
    Y = np.array( [1, 1, 1, 1, 1, 1.5, 1.5, 2, 2, 2, 2, 2, 1.5, 1.5] )
    
    A = np.transpose( np.vstack((
        eigenfunction_series(X, Y, n_terms, u11, mult),
        eigenfunction_series(X, Y, n_terms, u12, mult),
        eigenfunction_series(X, Y, n_terms, u13, mult),
        eigenfunction_series(X, Y, n_terms, u21, mult),
        eigenfunction_series(X, Y, n_terms, u22, mult),
        eigenfunction_series(X, Y, n_terms, u23, mult),
        eigenfunction_series(X, Y, n_terms, u31, mult),
        eigenfunction_series(X, Y, n_terms, u32, mult),
        eigenfunction_series(X, Y, n_terms, u33, mult),
        eigenfunction_series(X, Y, n_terms, u41, mult),
        eigenfunction_series(X, Y, n_terms, u42, mult),
        eigenfunction_series(X, Y, n_terms, u43, mult),
        eigenfunction_series(X, Y, n_terms, G( 1.5, 1.5 ), mult),
        eigenfunction_series(X, Y, n_terms, G( 7.5, 1.5 ), mult)
    )) )
    
    print(A)
    
    c = np.linalg.solve(A, np.array([1]*12 + [2]*2) )
    
    evaluate_fit_solution( c, Lx, Ly, n_terms )
    
    return c

c = fit_to_constant_along_lines()  
print( "Coefficients: ", c)  