import sys
from pathlib import Path
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))
from ConnectQuadToTet import save_to_obj
from surfacemesh import SurfaceMesh

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

def convert_boundary_txt_to_obj(obj_input_path, txt_path, obj_output_path, boundary = "Zero"):
    """
    Converts a txt file of vertex indices from an obj to a new obj file containing only the vertices and faces that are part of the specified boundary.
    obj_input_path: str - Path to the input OBJ file.
    txt_path: str - Path to the text file containing vertex indices.
    obj_output_path: str - Path where the output OBJ file will be saved.
    boundary: str - Specifies which boundary to filter by. Options are "Zero", "One", or "Both", default is "Zero".
    """
    def extract_indices(lines, label):
        """
        Extracts indices from lines of text based on a label.
        Splits lines by comma, grabbing every index on every line after the label until the next section header.
        lines: List[str] - Lines of text to search through.
        label: str - The label to look for in the lines.
        """
        collecting = False
        indices = []
        for line in lines:
            if line.startswith(f"# {label}:"):
                collecting = True
                line = line.split(':', 1)[1]  # Grab what's after the colon
            elif collecting and line.startswith('#'):
                break  # Stop at next section header
            if collecting:
                parts = [x.strip() for x in line.strip().split(',') if x.strip()]
                indices.extend(map(int, parts))
        return indices

    # Read the text file
    with open(txt_path, 'r') as f:
        lines = f.readlines()
        zero_boundary = extract_indices(lines, "Zero boundary verts")
        one_boundary = extract_indices(lines, "One boundary verts")
        if boundary == "Zero":
            valid_indices = set(zero_boundary)
        elif boundary == "One": 
            valid_indices = set(one_boundary)
        elif boundary == "Both":
            # Combine both sets of indices
            valid_indices = set(zero_boundary + one_boundary)
        else:
            valid_indices = set(zero_boundary)

    # Parse the .obj file
    vertices = []
    faces = []
    with open(obj_input_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                vertices.append(line.strip())
            elif line.startswith('f '):
                face = [int(part.split('/')[0]) for part in line.strip().split()[1:]]
                faces.append(face)

    # Filter faces with all vertices in valid set
    filtered_faces = [face for face in faces if all(idx in valid_indices for idx in face)]

    # Remap used vertex indices
    used_vertex_indices = sorted({idx for face in filtered_faces for idx in face})
    index_remap = {old: new for new, old in enumerate(used_vertex_indices, start=1)}

    # Get new vertex list
    new_vertices = [vertices[i - 1] for i in used_vertex_indices]  # 1-based indexing

    # Remap faces
    new_faces = [
        'f ' + ' '.join(str(index_remap[idx]) for idx in face)
        for face in filtered_faces
    ]

    # Write to new .obj file
    with open(obj_output_path, 'w') as f:
        for v in new_vertices:
            f.write(v + '\n')
        for face in new_faces:
            f.write(face + '\n')

    print(f"Filtered OBJ written to: {obj_output_path}")
