# A python script meant to run the process of taking a tet mesh and creating a quad mesh of one of the boundaries.
# The resulting vtk file path can be returned, along with a list of points and quads from the vtk file.

import sys
import subprocess
from pathlib import Path
import numpy as np
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
SWEEPS_DIR = SCRIPT_DIR.parent

path_to_api = Path(__file__).parent.parent / "build" / "src" / "api"
sys.path.insert(0, str(path_to_api))
import sweeps

def runQuadriflowProgram(input_file: str):
    """
    Runs the quadriflow program to create a quad mesh from a tet mesh input file.
    inputFile: str - Path to the tet mesh file.
    """
    quadriflow_executable = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "quadriflow"
    config_path = SWEEPS_DIR / "deps" / "ScaleUntrim" / "setting.config"
    if not quadriflow_executable.exists():
        raise FileNotFoundError(f"C++ executable not found at {quadriflow_executable}")
    
    command = [str(quadriflow_executable), input_file, config_path]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("OUTPUT:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"OUTPUT: {e.output}")
        print("END OF OUTPUT MESSAGE")
        print(f"ERROR OUTPUT: {e.stderr}")
        print("END OF ERROR MESSAGE")
        raise Exception("Error while running Quadriflow.")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception("Could not run Quadriflow.")

def generateQuadMesh(tri_mesh: sweeps.triMesh):
    # 1. turn objec tinto obj
    # Save the hook_base tri mesh to an OBJ file using Python's file io library
    tri_mesh_base_path = SCRIPT_DIR = Path(__file__).resolve().parent / "TetToQuadConnection" / "Input" / "MeshBase_new.obj"
    with open(tri_mesh_base_path, "w") as file:
        for vertex in mesh_base_tri.points:
            file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        for face in mesh_base_tri.tris:
            file.write(f"f {' '.join(str(idx + 1) for idx in face)}\n")

    # 2. turn obj into quad vtk
    # 3. make vtk an object



def generateQuadMeshOLD(tri_mesh: str):
    """
    A python function meant to take a tet mesh and create a quad mesh of the input as a vtk file. 
    Outputs the filepath of the new quad mesh.
    tetMesh: str - Path to the tet mesh file.
    """
    runQuadriflowProgram(tri_mesh)
    vtk_path = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "tempdir" / "quad.vtk"
    if not vtk_path.exists():
        raise FileNotFoundError("The generated quadrilateral vtk file does not exist.")
    return vtk_path

def createQuadMeshObjectFromVtk(vtk_path: str):
    """
    Takes a file path to a quad mesh as a vtk file and returns a sweeps.QuadMesh object derived from that vtk file.
    vtkPath: str - path to the vtk file to use in creating the sweeps.QuadMesh object.
    """
    f = open(vtk_path, "r")
    lines = f.readlines()
    f.close()
    num_points = -1
    points_start_index = -1
    num_quads = -1
    quads_start_index = -1
    for line in range(len(lines)):
        if lines[line].startswith("POINTS"):
            points_line = lines[line].split()
            num_points = int(points_line[1])
            points_start_index = line
        elif lines[line].startswith("CELLS"):
            quads_line = lines[line].split()
            num_quads = int(quads_line[1])
            quads_start_index = line
    if (num_points < 0) or (num_quads < 0) or (points_start_index < 0) or (quads_start_index < 0):
        raise Exception("Given Vtk file is empty or formatted incorrectly.")
    point_lines = lines[points_start_index+1:points_start_index+num_points+1]
    quad_lines = lines[quads_start_index+1:quads_start_index+num_quads+1]
    points = getPoints(point_lines)
    quads = getQuads(quad_lines)
    return sweeps.QuadMesh(points, quads)

def getPoints(lines):
    """
    Given lines of a file which contain vertices in vtk format, returns all vertices as a list of np arrays.
    lines: list - a list of lines, each line containing the 3 coordinates of a vertex.
    return: [np.array]
    """
    points = []
    for line in lines:
        nums = line.split()
        points.append(np.array(nums, dtype=float))
    return points

def getQuads(lines):
    """
    Given lines of a file which contain quadrilaterals in vtk format, returns all quads as a list of lists.
    lines: list = a list of lines, each line containing the number of vertices, and the 4 indices of the vertices in the quad.
    return: [list] - 2d list
    """
    quads = []
    for line in lines:
        nums = line.split()
        quads.append(list(map(int, nums[1:])))
    return quads

