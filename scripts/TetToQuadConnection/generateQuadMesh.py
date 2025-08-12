# A python script meant to run the process of taking a tet mesh and creating a quad mesh of one of the boundaries.
# It also outputs a file containing the barycentric coordinates of the quad vertices in the tet mesh faces.
import sys
import subprocess
from pathlib import Path
import os
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from src.ConnectQuadToTet import *
from src.PrepConnections import *

SCRIPT_DIR = Path(__file__).resolve().parent
SWEEPS_DIR = SCRIPT_DIR.parent.parent

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
    A python function meant to take a tet mesh and create a quad mesh of the input. 
    A obj and vtk are generated.
    Outputs the filepath of the new quad mesh.
    tetMesh: str - Path to the tet mesh file.
    """
    runQuadriflowProgram(tetMesh)
    vtkPath = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "tempdir" / "quad.vtk"
    if not vtkPath.exists():
        raise FileNotFoundError("The generated quadrilateral vtk file does not exist.")

    # Convert the vtk file to an obj file.
    quadMeshObjPath = SCRIPT_DIR / "output" / "quadMeshBoundary.obj"
    vtkToObj(vtkPath, quadMeshObjPath)
    if not quadMeshObjPath.exists():
        raise FileNotFoundError("The converted quadrilateral obj file does not exist.")

    # print results
    print("Quad Mesh of the boundary as an obj file is stored at: input/quadMeshBoundary.obj")
    print("The vtk version is stored at: sweeps/deps/ScaleUntrim/build/tempdir/quad.vtk")
