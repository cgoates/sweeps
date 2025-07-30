import sys
import argparse
from ConnectQuadToTet import *
from PrepConnections import *

# def main(tetMesh, quadMesh=None):
#     """
#     Take in a tet mesh, tet mesh start surface, and an optional quad mesh,
#     If no quad mesh is provided, create one with quadriflow.
#     create the quad surface using boundary and breadth first search.
#     superimpose the quad surface onto the tet surface.
#     """
#     return tetMesh, quadMesh

if __name__ == "__main__":
    """
    "Usage: python main.py <tetMesh> <quadMesh>"
    """
    parser = argparse.ArgumentParser(description="Connect a quad mesh to a tet mesh.")

    parser.add_argument("tetMesh", type=str, help="Path to the tet mesh file.")
    parser.add_argument("quadMesh", nargs='?', default=None, type=str, help="Path to the quad mesh file. If not provided, a quad mesh will be generated.")

    args = parser.parse_args()

    tetMesh = args.tetMesh
    quadMesh = args.quadMesh
    

    # add verification
    if tetMesh == None:
        print("Error: No tet mesh provided.")
        print("Usage: python main.py <tetMesh> <quadMesh>")
        print("If no quad mesh is provided, a quad mesh will be generated.")
        sys.exit(1)

    #This section is written in c++, so something has to be done to run it, unless I can find a way that works.
    #if quadriflow needs to be run, do it here.
    # if quadMesh == None:
    #     # run scaleUntrim quadriflow


    #     # convert the vtk to an obj

    #Now either way we have a tet and quad mesh.

    #Run the connection code.
    vtkToObj("../input/quad.vtk")

    connectQuadVertexOntoTriangle(tetMesh, quadMesh)

    #Return the desired output/ output the result
    #For now, both obj files, and the files containing which face each vertex resides in on the other mesh.
    #The files should be written, return the names of any created files.
    print("Output stored in output/quadVertices.txt")
