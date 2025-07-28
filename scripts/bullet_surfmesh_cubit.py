#!python
cubit.cmd('reset')
cubit.cmd('import "bullet_in_sphere.step"')
cubit.cmd('surf all except 1 copy')
cubit.cmd('delete vol 1')
cubit.cmd('unite surf all')
cubit.cmd('compress')
cubit.cmd('surf all size 2')
cubit.cmd('surf 4 scheme circle')
cubit.cmd('mesh surf all')

# Open output file for writing
with open('/Users/caleb/output.obj', 'w') as f:
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