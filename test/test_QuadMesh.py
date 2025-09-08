import sys
from pathlib import Path
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
PATH_TO_API = PROJECT_DIR / "build" / "src" / "api"
TEST_HOOK_FILE = SCRIPT_DIR / "data" / "hook.inp"
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PATH_TO_API))
from scripts.quadMeshingFunctions import generateQuadMesh
from scripts.scaleUntrimBuilder import buildScaleUntrim
import sweeps

def test_QuadMesh():
    """
    Run the process of Generating a quad mesh from the base of the hook mesh.
    The test will fail if any errors occur in the process.
    """
    try:
        # buildScaleUntrim()
        mesh = sweeps.loadFromFile(str(TEST_HOOK_FILE), "Surface12", "Surface10")
        tri_mesh_base = sweeps.baseOfSweep(mesh)
        quad_mesh = generateQuadMesh(tri_mesh_base)
        SWEEPS_DIR = Path(__file__).resolve().parent.parent
        vtk_path = SWEEPS_DIR / "deps" / "ScaleUntrim" / "build" / "tempDir" / "quad.vtk"
        print(SWEEPS_DIR)
        print(vtk_path)
        if not vtk_path.exists():
            raise FileNotFoundError("The generated quadrilateral vtk file does not exist.")
        u_values = np.linspace(0.0, 1.0, 30)
        mesh = sweeps.fitHexMeshToSweep(mesh, quad_mesh, u_values, debug=True)
    except Exception as e:
        raise AssertionError(f"Quad mesh couldn't be generated: {e}") from e

test_QuadMesh()