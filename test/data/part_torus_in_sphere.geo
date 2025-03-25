// Gmsh project created on Tue Mar 25 09:04:50 2025
//+
SetFactory("OpenCASCADE");
Torus(1) = {0, 0, 0, 0.5, 0.2, 6};
//+
Sphere(2) = {0, 0, 0, 1.2, -Pi/2, Pi/2, 2*Pi};

BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };
Mesh.MeshSizeFactor = 0.125;
Mesh 3;
