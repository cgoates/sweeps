from math import sin, cos, pi
import numpy as np
import sys
from pathlib import Path

path_to_api = Path(__file__).parent.parent.parent / "build" / "src" / "api"
print( path_to_api )
sys.path.insert(0, "/Users/caleb/sweeps/build/src/api")
import sweeps


def generate_points_in_circle(n, seed=42):
    """
    Generates n randomly distributed points within the unit circle
    
    Parameters:
    n (int): The number of points to generate
    seed (int): The random seed to use for reproducibility
    
    Returns:
    numpy.ndarray: An (n, 2) array of (x, y) coordinates of the points
    """
    np.random.seed(seed)

    # Generate random angles
    angles = np.random.uniform(0, 2 * np.pi, n)

    # Generate random distances from the origin
    distances = np.sqrt(np.random.uniform(0, 1, n))

    # Calculate (x, y) coordinates
    x = distances * np.cos(angles)
    y = distances * np.sin(angles)

    return [np.array(pt) for pt in zip(x, y)]


def parameterizeHook():
    # Load a tet mesh and the source and target surfaces from file. Surface12 is the source, and Surface10 is the target.
    hook = sweeps.loadFromFile(
        "/Users/caleb/sweeps/attempt-sweep/test/data/hook.inp", "Surface12", "Surface10")

    # This is a set of values at which you want to create level sets for the tracing.
    level_set_values = np.linspace(0.0, 1.0, 30)

    # this is a set of points in the unit circle that we will trace.
    trace_points = []
    for radius in [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        for theta in [0, 45, 90, 135, 180, 225, 270, 315]:
            trace_points.append(
                np.array([radius * cos(pi * theta / 180), radius * sin(pi * theta / 180)]))

    # perform the tracing and write out to files hook_level_sets.vtu and hook_traces.vtu which can be opened with paraview.
    sweeps.writeParameterizationToFile(
        hook, level_set_values, trace_points, "hook")


def parameterizeBunny():
    level_set_values = np.concatenate((np.linspace(0, 0.177, 35),
                                      np.linspace(0.177, 0.17904, 40),
                                      np.linspace(0.17904, 0.17908, 40),
                                      np.linspace(0.17908, 0.1796, 10),
                                      np.linspace(0.1796, 1.0, 20))).tolist()

    # 300 random points in the circle that we will trace.
    trace_points = generate_points_in_circle(300)

    bunny = sweeps.loadFromFile(
        "/Users/caleb/sweeps/attempt-sweep/test/data/stanford_bunny.inp", "placeholder", "placeholder")

    # The bunny doesn't have source and target sets defined in the inp file, so we create our own here.
    # The source is a set of vertices on the surface we are tracing from, and the target is a set of vertices on the surface we are tracing to.
    
    new_source = set()
    new_target = set()

    for [idx, pt] in enumerate(bunny.mesh.points):
        if pt[0] > 0.053:
            print( "Found source" )
            new_source.add(idx)
        elif pt[2] < -0.055:
            print( "Found Target" )
            new_target.add(idx)
            
    # bunny.source and bunny.target cannot be added to, they must be completely reassigned to be changed.
    # this is just a result of connecting python code to C++ code.
    bunny.source = new_source
    bunny.target = new_target
            
    print( "Source: ", bunny.source )
    print( "Target: ", bunny.target )

    sweeps.writeParameterizationToFile(
        bunny, level_set_values, trace_points, "bunny")


# Uncomment this to see what objects and functions are available from the sweeps module
# help( sweeps )

# Comment/uncomment these to run the different examples.
parameterizeHook()
# parameterizeBunny()
