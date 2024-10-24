edgeMesh = 10;
pointMesh = 1;

Point(1) = {-100, 100, 0, edgeMesh};
Point(2) = {-100, -100, 0, edgeMesh};
Point(3) = {100, -100, 0, edgeMesh};
Point(4) = {100, 100, 0, edgeMesh};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Point(11) = {-100, 100, -100, edgeMesh};
Point(12) = {-100, -100, -100, edgeMesh};
Point(13) = {100, -100, -100, edgeMesh};
Point(14) = {100, 100, -100, edgeMesh};
Line(11) = {12, 13};
Line(12) = {13, 14};
Line(13) = {14, 11};
Line(14) = {11, 12};
Curve Loop(11) = {11, 12, 13, 14};
Plane Surface(11) = {-11};
Line(15) = {12, 2};
Line(16) = {11, 1};
Line(17) = {14, 4};
Line(18) = {13, 3};
Curve Loop(12) = {15, 1, -18, -11};
Plane Surface(12) = {-12};
Curve Loop(13) = {12, 17, -2, -18};
Plane Surface(13) = {13};
Curve Loop(14) = {15, -4, -16, 14};
Plane Surface(14) = {14};
Curve Loop(15) = {13, 16, -3, -17};
Plane Surface(15) = {15};

Point(21) = {-3,-3, 0, pointMesh};
Point{21} In Surface{1};
Point(22) = {0,-3, 0, pointMesh};
Point{22} In Surface{1};
Point(23) = {3,-3, 0, pointMesh};
Point{23} In Surface{1};

Point(24) = {-3,0, 0, pointMesh};
Point{24} In Surface{1};
Point(25) = {0,0, 0, pointMesh};
Point{25} In Surface{1};
Point(26) = {3,0, 0, pointMesh};
Point{26} In Surface{1};

Point(27) = {-3,3, 0, pointMesh};
Point{27} In Surface{1};
Point(28) = {0,3, 0, pointMesh};
Point{28} In Surface{1};
Point(29) = {3,3, 0, pointMesh};
Point{29} In Surface{1};
Surface Loop(1) = {14, 12, 11, 15, 1, 13};
Volume(1) = {1};

Physical Surface("faces", 19) = {13, 11, 12, 14, 15};
Physical Volume("core", 20) = {1};
