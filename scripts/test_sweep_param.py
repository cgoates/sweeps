from math import sin, cos, pi
import numpy as np
import sys
from pathlib import Path

path_to_src = Path(__file__).parent.parent
path_to_api = path_to_src.parent / "build" / "src" / "api"
sys.path.insert(0, str(path_to_api))
import sweeps

def meshHook(single_patch=True):
    # Load a tet mesh and the source and target surfaces from file. Surface12 is the source, and Surface10 is the target.
    hook = sweeps.loadFromFile(
        str( path_to_src / "test" / "data" / "hook.inp" ), "Surface12", "Surface10")

    # This is a set of u values at which you want the mesh to have points.
    # Try changing this to get an idea for what it does.  The values should all be between 0 and 1.
    u_values = np.linspace(0.0, 1.0, 30)

    # for single patch, the resulting mesh will have n_elems_st x n_elems_st x len( u_values ) elements.
    # for the five patch case, it will be 5 times that amount.
    n_elems_st = 5

    # Setting the debug flag to True will check the mesh jacobians, and write out a vtk file to visualize the mesh and any elements with negative jacobians.
    mesh = sweeps.fitSinglePatchHexMeshToSweep(hook, n_elems_st, u_values, debug=True) if single_patch else sweeps.fitFivePatchHexMeshToSweep(hook, n_elems_st, u_values, debug=True)

    # The resulting hex mesh has a list of points and a list of hexes.
    # The API documentation has this description of the data structure:
    # A simple hex mesh, with a list of points, and a list of hexes.  The vertices of the hexes are
    # ordered by their coordinates in the hex-local coordinates system as (0,0,0), (1,0,0), (0,1,0),
    # (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1).
    print( len( mesh.points ) )
    print( len( mesh.hexes ) )
    print( [ mesh.points[i] for i in mesh.hexes[0] ] )

meshHook(single_patch=False)
