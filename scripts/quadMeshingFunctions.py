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

def _writeConfigWithMagnitude(magnitude_factor: float, angle_for_sharp: float = None) -> Path:
    """
    Write a temporary config file with the given magnitude_factor (and optionally
    angle_for_sharp) and return its path.
    """
    import tempfile
    base_config = SWEEPS_DIR / "deps" / "ScaleUntrim" / "setting.config"
    tmp = Path(tempfile.mktemp(suffix=".config"))
    lines = base_config.read_text().splitlines()
    with open(tmp, "w") as f:
        for line in lines:
            if line.startswith("magnitude_factor:"):
                f.write(f"magnitude_factor: {magnitude_factor}\n")
            elif angle_for_sharp is not None and line.startswith("angle_for_sharp:"):
                f.write(f"angle_for_sharp: {angle_for_sharp}\n")
            else:
                f.write(line + "\n")
    return tmp

def runQuadriflowProgram(input_file: str, config_path: Path = None, sharp: bool = False):
    """
    Runs the quadriflow program to create a quad mesh from a tet mesh input file.
    input_file: str - Path to the tet mesh file.
    config_path: Path - Optional path to config file. Uses default setting.config if None.
    sharp: bool - If True, passes --sharp to enable sharp feature preservation.
    """
    quadriflow_executable = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "QuadriFlow"
    if config_path is None:
        config_path = SWEEPS_DIR / "deps" / "ScaleUntrim" / "setting.config"
    if not quadriflow_executable.exists():
        raise FileNotFoundError(f"C++ executable not found at {quadriflow_executable}")

    command = [str(quadriflow_executable), str(input_file), str(config_path),
               "--sharp" if sharp else "--no-sharp"]
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

def generateQuadMesh(tri_mesh: sweeps.TriMesh, sharp_angle: float = 0):
    """
    Takes a sweeps.TriMesh object and creates an obj using the data. The obj is used to create a quad mesh vtk file.
    The vtk file is used to generate a sweeps.QuadMesh object which is then returned.
    If QuadriFlow encounters numerical issues (NaN), automatically retries with progressively larger
    magnitude_factor values (coarser quad mesh) until it succeeds.
    tri_mesh: sweeps.TriMesh - A TriMesh object containing points and face data for the mesh to be generated.
    sharp_angle: float - If > 0, enables sharp feature preservation using this value as the dihedral
                         angle threshold (in degrees). Edges where adjacent faces differ by more than
                         this angle are treated as sharp features. If <= 0 (default), sharp feature
                         preservation is disabled.
    Returns: a sweeps.QuadMesh object containing the points and face (quads) data from the vtk output.
    """
    tri_mesh_base_path = Path(__file__).resolve().parent / "MeshBase_new.obj"
    saveTriMeshAsObj(tri_mesh, str(tri_mesh_base_path))

    vtk_path = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "tempDir" / "quad.vtk"
    use_sharp = sharp_angle > 0

    # Try progressively larger magnitude_factor values until QuadriFlow succeeds.
    # A higher magnitude_factor produces a coarser quad mesh but avoids NaN failures
    # that can occur with fine meshes due to degenerate configurations in the solver.
    magnitude_factors = [0.89,1,1.8, 2, 2.3, 3, 5, 6, 8, 10]
    for mf in magnitude_factors:
        if vtk_path.exists():
            vtk_path.unlink()
        tmp_config = _writeConfigWithMagnitude(mf, angle_for_sharp=sharp_angle if use_sharp else None)
        try:
            print(f"Running QuadriFlow with magnitude_factor={mf}, sharp={'on (angle=' + str(sharp_angle) + ')' if use_sharp else 'off'}...")
            runQuadriflowProgram(tri_mesh_base_path, tmp_config, sharp=use_sharp)
        finally:
            tmp_config.unlink(missing_ok=True)
        if vtk_path.exists():
            print(f"QuadriFlow succeeded with magnitude_factor={mf}.")
            break
        print(f"QuadriFlow did not produce output with magnitude_factor={mf}, retrying...")
    else:
        raise RuntimeError(
            f"QuadriFlow failed to produce a quad mesh after trying magnitude_factors={magnitude_factors}. "
            "The input triangle mesh may be too complex or degenerate for this solver."
        )

    return createQuadMeshObjectFromVtk(str(vtk_path))

def createQuadMeshObjectFromVtk(vtk_path: str):
    """
    Takes a file path to a quad mesh as a vtk file and returns a sweeps.QuadMesh object derived from that vtk file.
    vtk_path: str - path to the vtk file to use in creating the sweeps.QuadMesh object.
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

def saveTetMeshAsVtu(tet_mesh: sweeps.TetMesh, path_to_file: str):
    """
    Writes a sweeps.TetMesh to a VTU file (VTK XML unstructured grid).
    tet_mesh: sweeps.TetMesh - a TetMesh object with .points and .tets.
    path_to_file: str - destination .vtu file path.
    """
    points = tet_mesh.points
    tets = [tet.vertices for tet in tet_mesh.tets]
    n_pts = len(points)
    n_cells = len(tets)
    with open(path_to_file, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n')
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for p in points:
            f.write(f'          {p[0]} {p[1]} {p[2]}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        for tet in tets:
            f.write(f'          {" ".join(str(v) for v in tet)}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        for i in range(1, n_cells + 1):
            f.write(f'          {i * 4}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in tets:
            f.write('          10\n')  # VTK_TETRA = 10
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')


def saveQuadMeshAsVtu(quad_mesh: sweeps.QuadMesh, path_to_file: str):
    """
    Writes a sweeps.QuadMesh to a VTU file (VTK XML unstructured grid).
    quad_mesh: sweeps.QuadMesh - a QuadMesh object with .points and .quads.
    path_to_file: str - destination .vtu file path.
    """
    points = quad_mesh.points
    quads = quad_mesh.quads
    n_pts = len(points)
    n_cells = len(quads)
    with open(path_to_file, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n')
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for p in points:
            f.write(f'          {p[0]} {p[1]} {p[2]}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        for quad in quads:
            f.write(f'          {" ".join(str(v) for v in quad)}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        for i in range(1, n_cells + 1):
            f.write(f'          {i * 4}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in quads:
            f.write('          9\n')  # VTK_QUAD = 9
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')
        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')


def saveTriMeshAsObj(mesh_base_tri: sweeps.TriMesh, path_to_file: str):
    """
    Takes in a sweeps.TriMesh object and a desired file path to save to and
    creates an obj file using the data from the TriMesh object.
    mesh_base_tri: sweeps.TriMesh - a TriMesh object created using the sweeps api.
    path_to_file: str - the path for the obj file to be saved to.
    """
    with open(path_to_file, "w") as file:
        for vertex in mesh_base_tri.points:
            file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        for face in mesh_base_tri.tris:
            file.write(f"f {' '.join(str(idx + 1) for idx in face)}\n")
