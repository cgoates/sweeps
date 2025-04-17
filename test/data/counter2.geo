SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, -0.5, 0, 0, 1, 1.5, 2*Pi};
Cylinder(2) = {1, -0.5, 0, 0, 1, 0, 1.5, 2*Pi};
BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }
Mesh.MeshSizeFactor = 0.15;
Mesh 3;
