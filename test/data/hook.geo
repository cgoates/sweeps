SetFactory("OpenCASCADE");
v() = ShapeFromFile("hook.iges");
Mesh.MeshSizeFactor = 0.25/3;
Mesh 3;//+
//Physical Surface("source", 7) = {3};
