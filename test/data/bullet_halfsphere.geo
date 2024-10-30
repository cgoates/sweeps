SetFactory("OpenCASCADE");
v() = ShapeFromFile("bullet_halfsphere.stp");
bbox() = BoundingBox Volume{v()};
xmin = bbox(0);
ymin = bbox(1);
zmin = bbox(2);
xmax = bbox(3);
ymax = bbox(4);
zmax = bbox(5);
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };
//Box(50)= {0,-40,-40,40,40,40};
//Translate{-20,0,20}{Volume{50};};
//BooleanDifference{ Volume{v()}; Delete; }{ Volume{50}; Delete; };
Mesh.MeshSizeFactor = 0.25;
//Mesh 3;//+
//Physical Surface("source", 7) = {3};
