// Dodecahedron generation in Gmsh

// Parameters
lc = 1.5; // Characteristic length for mesh size

// Golden ratio
phi = (1 + Sqrt(5)) / 2;

// Define points
Point(1) = {-1, 1, 1, lc};
Point(2) = {1, 1, 1, lc};
Point(3) = {1, -1, 1, lc};
Point(4) = {-1, -1, 1, lc};
Point(5) = {-1, 1, -1, lc};
Point(6) = {1, 1, -1, lc};
Point(7) = {1, -1, -1, lc};
Point(8) = {-1, -1, -1, lc};
Point(9) = {0, phi, 1/phi, lc};
Point(10) = {0, phi, -1/phi, lc};
Point(11) = {0, -phi, 1/phi, lc};
Point(12) = {0, -phi, -1/phi, lc};
Point(13) = {1/phi, 0, phi, lc};
Point(14) = {-1/phi, 0, phi, lc};
Point(15) = {1/phi, 0, -phi, lc};
Point(16) = {-1/phi, 0, -phi, lc};
Point(17) = {phi, 1/phi, 0, lc};
Point(18) = {phi, -1/phi, 0, lc};
Point(19) = {-phi, 1/phi, 0, lc};
Point(20) = {-phi, -1/phi, 0, lc};

// Define lines
Line(1) = {1, 9};   Line(2) = {9, 2};   Line(3) = {2, 13};  Line(4) = {13, 3};  Line(5) = {3, 11};
Line(6) = {11, 4};  Line(7) = {4, 14};  Line(8) = {14, 1};  Line(9) = {1, 19};  Line(10) = {19, 5};
Line(11) = {5, 10}; Line(12) = {10, 9}; Line(13) = {2, 17}; Line(14) = {17, 6}; Line(15) = {6, 10};
Line(16) = {5, 16}; Line(17) = {16, 8}; Line(18) = {8, 20}; Line(19) = {20, 19};Line(20) = {4, 20};
Line(21) = {3, 18}; Line(22) = {18, 7}; Line(23) = {7, 12}; Line(24) = {12, 11};Line(25) = {6, 15};
Line(26) = {15, 7}; Line(27) = {8, 12}; Line(28) = {16, 15};Line(29) = {17, 18};Line(30) = {13, 14};

// // Define surfaces (faces)
Line Loop(1) = {1, 2, 3, 30, 8};   Plane Surface(1) = {1};
Line Loop(2) = {9, 10, 11, 12, -1};        Plane Surface(2) = {2};
Line Loop(3) = {2, 13, 14, 15, 12};       Plane Surface(3) = {3};
Line Loop(4) = {3, 4, 21, -29, -13};       Plane Surface(4) = {4};
Line Loop(5) = {7, 8, 9, -19, -20}; Plane Surface(5) = {5};
Line Loop(6) = {4, 5, 6, 7, -30};      Plane Surface(6) = {6};
Line Loop(7) = {6, 20, -18, 27, 24};           Plane Surface(7) = {7};
Line Loop(8) = {10, 16, 17, 18, 19};         Plane Surface(8) = {8};
Line Loop(9) = {23, 24, -5, 21, 22};      Plane Surface(9) = {9};
Line Loop(10) = {25, 26, -22, -29, 14};    Plane Surface(10) = {10};
Line Loop(11) = {28, -25, 15, -11, 16};   Plane Surface(11) = {11};
Line Loop(12) = {28, 26, 23, -27, -17};   Plane Surface(12) = {12};

// Create surface loop and volume
Surface Loop(1) = {1:12};
Volume(1) = {1};

// Set mesh size
Mesh.MeshSizeFactor = lc;


// Generate 3D mesh
Mesh 3;