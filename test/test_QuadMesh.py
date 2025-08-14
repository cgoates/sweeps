import sys
from pathlib import Path
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
PATH_TO_API = PROJECT_DIR / "build" / "src" / "api"
TEST_HOOK_FILE = SCRIPT_DIR / "data" / "hook.inp"
TEST_HOOK_BASE_FILE = SCRIPT_DIR / "data" / "hookBaseTri.obj"
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PATH_TO_API))
from scripts.TetToQuadConnection.quadMeshingFunctions import generateQuadMesh, createQuadMeshObjectFromVTK
import sweeps

def test_QuadMesh():
    """
    Run the process of Generating a quad mesh from the base of the hook mesh.
    The test will fail if any errors occur in the process.
    """
    try:
        mesh = sweeps.loadFromFile(str(TEST_HOOK_FILE), "Surface12", "Surface10")
        mesh_path = generateQuadMesh(tetMesh=str(TEST_HOOK_BASE_FILE))
        quad_mesh = createQuadMeshObjectFromVTK(mesh_path)
        u_values = np.linspace(0.0, 1.0, 30)
        mesh = sweeps.fitHexMeshToSweep(mesh, quad_mesh, u_values, debug=True)
    except:
        raise AssertionError("Quad mesh couldn't be generated.")