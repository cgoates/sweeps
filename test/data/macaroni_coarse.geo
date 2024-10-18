SetFactory("OpenCASCADE");
v() = ShapeFromFile("TorusPipe.stp");
Mesh.MeshSizeFactor = 0.25;
Mesh 3;//+
//Physical Surface("source", 7) = {3};
