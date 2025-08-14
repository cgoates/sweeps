# A python script meant to run the process of taking a tet mesh and creating a quad mesh of one of the boundaries.
# The resulting vtk file path can be returned, along with a list of points and quads from the vtk file.

import sys
import subprocess
from pathlib import Path
import numpy as np
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
SWEEPS_DIR = SCRIPT_DIR.parent.parent

path_to_api = Path(__file__).parent.parent.parent / "build" / "src" / "api"
sys.path.insert(0, str(path_to_api))
import sweeps

def runQuadriflowProgram(inputFile: str):
    """
    Runs the quadriflow program to create a quad mesh from a tet mesh input file.
    inputFile: str - Path to the tet mesh file.
    """
    quadriflowExecutable = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "quadriflow"
    configPath = SWEEPS_DIR / "deps" / "ScaleUntrim" / "setting.config"
    
    if not quadriflowExecutable.exists():
        raise FileNotFoundError(f"C++ executable not found at {quadriflowExecutable}")
    
    
    command = [str(quadriflowExecutable), inputFile, configPath]

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
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

def generateQuadMesh(tetMesh: str):
    """
    A python function meant to take a tet mesh and create a quad mesh of the input as a vtk file. 
    Outputs the filepath of the new quad mesh.
    tetMesh: str - Path to the tet mesh file.
    """
    runQuadriflowProgram(tetMesh)
    vtkPath = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "tempdir" / "quad.vtk"
    if not vtkPath.exists():
        raise FileNotFoundError("The generated quadrilateral vtk file does not exist.")
    return vtkPath

def createQuadMeshObjectFromVTK(vtkPath: str):
    """
    Takes a file path to a quad mesh as a vtk file and returns a sweeps.QuadMesh object derived from that vtk file.
    vtkPath: str - path to the vtk file to use in creating the sweeps.QuadMesh object.
    """
    f = open(vtkPath, "r")
    lines = f.readlines()
    f.close()
    numPoints = -1
    pointsStartIndex = -1
    numQuads = -1
    quadsStartIndex = -1
    for line in range(len(lines)):
        if lines[line].startswith("POINTS"):
            pointsLine = lines[line].split()
            numPoints = int(pointsLine[1])
            pointsStartIndex = line
        elif lines[line].startswith("CELLS"):
            quadsLine = lines[line].split()
            numQuads = int(quadsLine[1])
            quadsStartIndex = line
    if (numPoints < 0) or (numQuads < 0) or (pointsStartIndex < 0) or (quadsStartIndex < 0):
        raise Exception("Given Vtk file is empty or formatted incorrectly.")
    pointLines = lines[pointsStartIndex+1:pointsStartIndex+numPoints+1]
    quadLines = lines[quadsStartIndex+1:quadsStartIndex+numQuads+1]
    points = getPoints(pointLines)
    quads = getQuads(quadLines)
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