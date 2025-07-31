import math
from numpy import *
import builtins  # For explicit use of Python's built-in min and max
from pathlib import Path
import sys
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))
from surfacemesh import SurfaceMesh

def project_point_onto_triangle_plane(point, v0, v1, v2):
    """
    Projects a 3D point onto a triangular face defined by vertices v0, v1, v2.
    
    Args:
        point (np.array): The 3D point to project (x, y, z).
        v0, v1, v2 (np.array): The vertices of the triangle (x, y, z).
        
    Returns:
        np.array: The projected point on the plane of the triangle.
    """
    # Convert to numpy arrays
    point = array(point)
    v0, v1, v2 = array(v0), array(v1), array(v2)
    
    # Calculate plane normal
    normal = cross(v1 - v0, v2 - v0)
    normal = normal / linalg.norm(normal)
    
    # Find the projection of the point onto the plane
    point_to_v0 = point - v0
    distance_to_plane = dot(point_to_v0, normal)
    projected_point = point - distance_to_plane * normal

    return  projected_point

def project_point_onto_triangle_edge(point, v0, v1):
    """
    Projects a 3D point onto an edge defined by vertices v0, v1.
    
    Args:
        point (np.array): The 3D point to project (x, y, z).
        v0, v1 (np.array): The vertices of the edge (x, y, z).
        
    Returns:
        np.array: The projected point on the edge.
        float: The distance from the original point to the projected point.
    """
    # Convert to numpy arrays
    point = array(point)
    v0, v1 = array(v0), array(v1)
    
    # Calculate edge vector
    edge = v1 - v0
    edge_length = linalg.norm(edge)
    edge = edge / edge_length  # Normalize edge vector
    
    # Find the projection of the point onto the edge
    point_to_v0 = point - v0
    distance_to_edge = dot(point_to_v0, edge)

    distance_to_edge_float = float(distance_to_edge)
    edge_length_float = float(edge_length)

    # Clamp the projection distance to ensure it's within the edge segment
    distance_to_edge = builtins.max(0, builtins.min(distance_to_edge_float, edge_length_float)) # throwing type error

    projected_point = v0 + distance_to_edge * edge

    # Now we need to return the closest point on the edge itself, along with the distance from the point we started with.
    distance = linalg.norm(point - projected_point)

    return  projected_point , distance


def convertToBarycentric(p, a, b, c):
    """
    Compute barycentric coordinates (u, v, w) for
    point p with respect to triangle (a, b, c)
    Args:
        p (np.array): The point to convert (x, y, z).
        a, b, c (np.array): The vertices of the triangle (x, y, z).
    Returns:
        np.array: Barycentric coordinates (u, v, w) of point p with respect to triangle (a, b, c).
    """
    v0 = b - a
    v1 = c - a
    v2 = p - a 
    d00 = dot(v0, v0)
    d01 = dot(v0, v1)
    d11 = dot(v1, v1)
    d20 = dot(v2, v0)
    d21 = dot(v2, v1)
    denom = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1.0 - v - w
    return array([u, v, w])

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
        #obj_file.write("f 7 9 11 0\n")

def calcClosestVertice(triMesh, chosenVertex):
    """
    Finds the vertex on the triangular mesh with the smallest distance to the chosen vertex from the quad mesh.
    triMesh: SurfaceMesh - The triangle mesh object.
    chosenVertex: np.array - The vertex from the quadrilateral mesh to find the closest triangle vertex to.
    """
    minDistance = linalg.norm(chosenVertex - triMesh.vs[0])
    index = 0
    vertexIndex = 0
    for vertex in triMesh.vs:
            #calculate the distance between the chosen vertex and the current vertex
        distance = linalg.norm(chosenVertex - vertex)
            #track the vertex with min distance
        if distance < minDistance:
            minDistance = distance
            vertexIndex = index
        index += 1
    return vertexIndex

def findResidingFace(triMesh, chosenVertex, vertexIndicesOfAdjacentFaces):
    """
    Finds the face on the triangle mesh object which the chosen vertex resides in, if any.
    triMesh: SurfaceMesh - The triangle mesh object.
    chosenVertex: np.array - The vertex from the quadrilateral mesh to check against the triangle faces.
    vertexIndicesOfAdjacentFaces: 2D list - A list of lists containing the indices of vertices for each face adjacent to the closest triangle vertex of the chosen vertex.
    """
    badCoord = True
    closestFaceIndex = -1
    for face in vertexIndicesOfAdjacentFaces:
        closestFaceIndex += 1
        chosenTriFaceVertices = []
        for vert in face:
            chosenTriFaceVertices.append(triMesh.vs[vert])

        #Project to closest point on plane of a triangle.
        #project the quad point onto the triangle
        projectedPoint = project_point_onto_triangle_plane(chosenVertex, chosenTriFaceVertices[0], chosenTriFaceVertices[1], chosenTriFaceVertices[2])

        #use barycentric coordinates to express it.
        barycentric = convertToBarycentric(projectedPoint, chosenTriFaceVertices[0], chosenTriFaceVertices[1], chosenTriFaceVertices[2])
        if barycentric[0] >= 0 and barycentric[0] <= 1 and barycentric[1] >= 0 and barycentric[1] <= 1 and barycentric[2] >= 0 and barycentric[2] <= 1:
            badCoord = False
            break
    return badCoord, projectedPoint, barycentric, closestFaceIndex

def findClosestFace(triMesh, facesAdjacentToVertex, projectedPoint):
    """
    Finds the closest face on the triangle mesh to the projected point.
    triMesh: SurfaceMesh - The triangle mesh object.
    facesAdjacentToVertex: list - A list of indices of faces adjacent to the closest triangle vertex of the chosen vertex.
    projectedPoint: np.array - The point projected onto the triangle plane.
    """
    faceIndex = 0
    minDistance = math.inf
    for face in facesAdjacentToVertex:
                # get the halfedges of the face.
                # use to_vertex and opposite_edge -> to_vertex to get the vertices of the edges.
                # Do it for all edges until at the first. 
                # Keep the smallest distance then move to next face.
        edges = triMesh.get_face_halfedges(face)
        edgeIndex = 0
        for edge in edges:
            v0 = triMesh.vs[edge.ToVertex()]
            v1 = triMesh.vs[triMesh.halfedges[edge.next_he].ToVertex()]
            projectedPointOnEdge , distance = project_point_onto_triangle_edge(projectedPoint, v0, v1)
            if distance < minDistance:
                minDistance = distance
                closestFaceIndex = faceIndex
                closestProjectedPointOnEdge = projectedPointOnEdge
            edgeIndex += 1
        faceIndex += 1
    return closestProjectedPointOnEdge, closestFaceIndex

def connectQuadVertexOntoTriangle(triMeshFile, quadMeshFile, outputPath):
    """
    Takes a triangle mesh and a quadrilateral mesh, uses the SurfaceMesh library to make half-edge data structures for both, 
    finds the closest triangle face for each quadrilateral vertex,
    returns a txt file containing the barycentric coordinates of each quad mesh vertex with respect to the closest triangle face.
    triMeshFile: str - Path to the triangle mesh file (OBJ format).
    quadMeshFile: str - Path to the quadrilateral mesh file (OBJ format).
    outputPath: str - Path to the output text file where results will be saved.
    """
    # create a mesh of both files.
    triMesh = SurfaceMesh.FromOBJ_FileName(triMeshFile)
    quadMesh = SurfaceMesh.FromOBJ_FileName(quadMeshFile)

    quadFile = open(outputPath, "w")
    quadFile.write("Quadrilateral Vertex Index, Triangle Face Index, Barycentric Coordinates (u v w)\n")

    #determine how many quad verices there are.
    # for loop that many times 
    numQuadVertices = len(quadMesh.vs)
    for i in range(numQuadVertices):
        chosenVertex = quadMesh.vs[i]
        vertexIndex = calcClosestVertice(triMesh, chosenVertex)
        facesAdjacentToVertex = triMesh.GetFacesAdjacentToVertex(vertexIndex) # these are the indices of the faces adjacent to the closest triangle vertex.

        vertexIndicesOfAdjacentFaces = []
        for face in facesAdjacentToVertex:
            vertexIndicesOfAdjacentFaces.append(triMesh.faces[face]) # 2D array

        # If the point resides in a face when projected onto the face's plane, find out the barycentric coordinates 
        badCoord, projectedPoint, barycentric, closestFaceIndex = findResidingFace(triMesh, chosenVertex, vertexIndicesOfAdjacentFaces)
        #if none of them had barycentric coordinates all between 0 and 1, find the closest face and use its barycentric coordinates.
        if badCoord:
            closestProjectedPointOnEdge, closestFaceIndex = findClosestFace(triMesh, facesAdjacentToVertex, projectedPoint)
            chosenTriFaceVertices = []
            for vert in triMesh.faces[closestFaceIndex]:
                chosenTriFaceVertices.append(triMesh.vs[vert])
            barycentric = convertToBarycentric(closestProjectedPointOnEdge, chosenTriFaceVertices[0], chosenTriFaceVertices[1], chosenTriFaceVertices[2])
        
        #Store the triangle and barycentric coordinates.
        closestFaceMeshIndex = facesAdjacentToVertex[closestFaceIndex]
        #chosenTriFaceVertices - the triangle face vertices list
        #barycentric - the barycentric coordinates of the projected point

        #write the vertice each iteration of the for loop.
        # Format: "the current vertex indice from the quadrilateral mesh, the face indice it resides in or is closest to, the barycentric coordinates of the vertice in that face."
        quadFile.write(f"{i}, {closestFaceMeshIndex}, {barycentric[0]} {barycentric[1]} {barycentric[2]}\n")
    quadFile.close()
