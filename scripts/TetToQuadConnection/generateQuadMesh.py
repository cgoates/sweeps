# A python script meant to run the process of taking a tet mesh and creating a quad mesh of one of the boundaries.
# It also outputs a file containing the barycentric coordinates of the quad vertices in the tet mesh faces.
import sys
import subprocess
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from TetToQuadConnection.QuadMeshClass import QuadMesh


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
    # TODO: Uncomment the next line and create a quadMesh object, remove all lines following it
    # return QuadMesh(vtkPath)
    quadMeshObjPath = SCRIPT_DIR / "output" / "quadMeshBoundary.obj"
    vtkToObj(vtkPath, quadMeshObjPath)
    if not quadMeshObjPath.exists():
        raise FileNotFoundError("The converted quadrilateral obj file does not exist.")

####################
# TODO: Delete all code below this once the vtk to quadMesh object conversion is implemented, removing the need for obj files.
####################
def save_to_obj(file_path, vertices, faces):
    """
    Saves vertices and faces to a .obj file.

    Args:
        file_path (str): The output file path.
        vertices (np.array): A (N, 3) array of vertex positions (x, y, z).
        faces (np.array): A (M, 3) or (M, 4) array of faces (vertex indices, 1-based).
    """
    with open(file_path, 'w') as obj_file:
        # Write vertices
        for vertex in vertices:
            obj_file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        
        # Write faces
        for face in faces:
            face_str = " ".join(str(int(index) + 1) for index in face)  # Convert to 1-based indexing
            obj_file.write(f"f {face_str}\n")


def vtkToObj(vtkMesh, filename = "../input/convertedVTK.obj"):
    """
    Converts a VTK mesh file to an OBJ file format.
    vtkMesh: str - Path to the VTK mesh file.
    filename: str - Path where the converted OBJ file will be saved.
    """
    # get the vertices and faces from the vtk mesh.
    # obj is 1 indexed. vtk is 0 indexed.
    with open(vtkMesh, "r") as quadFile: 
        lines = quadFile.readlines()
        vertices = []
        faces = []
        numPoints = 0
        numCells = 0
        for line in lines:
            #startwith POINTS : vertices : x y z
            if numPoints > 0:
                vertices.append((line.split()))
                numPoints -= 1
            if line.startswith("POINTS"):
                numPoints = int(line.split()[1])     
            #startwith CELLS : Faces : 4 vertex vertex vertex vertex : the vertex is an index number, same as .obj
            if numCells > 0:
                nums = line.split()[1:]
                for i in range(len(nums)):
                    nums[i] = str(nums[i])
                faces.append(nums)
                numCells -= 1
            if line.startswith("CELLS"):
                numCells = int(line.split()[1])

    #write them to the obj
    save_to_obj(f"{filename}", vertices, faces)
    return None
