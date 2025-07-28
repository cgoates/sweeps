#!python
cubit.cmd('reset')
cubit.cmd('import step "APOLLO_capsule_immersed.step" heal')
cubit.cmd('vol 1 size 200')
cubit.cmd('mesh surf 2 3 4 5')
cubit.cmd('surf 2 3 4 5 copy')
cubit.cmd('delete vol 1')
cubit.cmd('compress')
cubit.cmd('webcut surf 2 3 with plane yplane')
cubit.cmd('imprint body all')
cubit.cmd('merge body all')
cubit.cmd('mesh surf all')
cubit.cmd('delete mesh surf 4')
cubit.cmd('copy mesh surface 1 onto  surface 4  source curve 12  source vertex 1  target curve 17  target vertex 14   smooth ')
cubit.cmd('surface 4  smooth scheme centroid area pull ')
cubit.cmd('smooth surface 4')
cubit.cmd('compress')

# Open output file for writing
with open('/Users/caleb/capsule.obj', 'w') as f:
    # Get and write vertices
    nodes = cubit.get_entities("node")
    for node in nodes:
        coords = cubit.get_nodal_coordinates(node)
        f.write(f"v {coords[0]} {coords[1]} {coords[2]}\n")
    
    f.write("\n")  # Blank line for readability
    
    # Get and write faces
    faces = cubit.get_entities("face")
    print(f"Processing {len(faces)} faces...")
    
    for face in faces:
        connectivity = cubit.get_connectivity("face", face)
        # OBJ face format: f v1 v2 v3 ... (1-indexed)
        face_line = "f " + " ".join(str(node_id) for node_id in connectivity)
        f.write(face_line + "\n")

print("OBJ file 'output.obj' has been created successfully!")