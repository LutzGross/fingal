// GMSH geometry file
// 
Mesh.MshFileVersion = 2.2;
// ... electrodes .....
numElectrodes = 32;
distanceElectrodes = 3;
lengthLine = (numElectrodes-1) * distanceElectrodes;
// ... core region ....
depthDomain  = lengthLine * 0.3;
widthY = lengthLine * 0.3;
widthX = lengthLine + widthY;

// ... Anomaly ...
anomalyDepth = 4 * distanceElectrodes;
anomalyLengthX = 10 * distanceElectrodes;
anomalyLengthY = 5 * distanceElectrodes;
anomalyLengthZ = 5 * distanceElectrodes;

// ... Padding ......
paddingWidth = 0.6 * lengthLine;
boundingBoxX =  widthX + 2 * paddingWidth;
boundingBoxY =  widthY + 2 * paddingWidth;
boundingBoxZ =  depthDomain + paddingWidth;


// ... element sizes ...
meshSizeEdges = distanceElectrodes;
meshSizeElectrodes = distanceElectrodes/4;
meshSizeBoundingBox = boundingBoxY/10;

// .... hexahedral anomaly ........................................................................
Point(9) = { 0, -anomalyLengthY/2, -anomalyDepth, 		meshSizeEdges};
Point(10) = { anomalyLengthX, -anomalyLengthY/2, -anomalyDepth, meshSizeEdges};
Point(11) = {0,  anomalyLengthY/2, -anomalyDepth, 		meshSizeEdges};
Point(12) = { anomalyLengthX,  anomalyLengthY/2, -anomalyDepth, meshSizeEdges};

Point(13) = {0, -anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 		meshSizeEdges};
Point(14) = { anomalyLengthX, -anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 	meshSizeEdges};
Point(15) = {0,  anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 		meshSizeEdges};
Point(16) = { anomalyLengthX,  anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 	meshSizeEdges};

Line(13) = {9, 10};
Line(14) = {10, 12};
Line(15) = {12, 11};
Line(16) = {11, 9};
Line(17) = {13, 14};
Line(18) = {13, 14};
Line(19) = {14, 16};
Line(20) = {16, 15};
Line(21) = {15, 13};
Line(22) = {13, 9};
Line(23) = {14, 10};
Line(24) = {16, 12};
Line(25) = {15, 11};
Curve Loop(7) = {17, 23, -13, -22};
Plane Surface(7) = {7};
Curve Loop(8) = {19, 24, -14, -23};
Plane Surface(8) = {8};
Curve Loop(9) = {20, 25, -15, -24};
Plane Surface(9) = {9};
Curve Loop(10) = {16, 13, 14, 15};
Plane Surface(10) = {10};
Curve Loop(11) = {25, 16, -22, -21};
Plane Surface(11) = {11};
Curve Loop(12) = {17, 19, 20, 21};
Plane Surface(12) = {12};
Surface Loop(3) = {11, 9, 12, 7, 8, 10};
Volume(2) = {3};
// ........ End Anomaly ...................

// ... Core .................................................

Point(1) = {-widthX/2, -widthY/2, 0, meshSizeEdges};
Point(2) = { widthX/2, -widthY/2, 0, meshSizeEdges};
Point(3) = {-widthX/2,  widthY/2, 0, meshSizeEdges};
Point(4) = { widthX/2,  widthY/2, 0, meshSizeEdges};

Point(5) = {-widthX/2, -widthY/2, -depthDomain, meshSizeEdges};
Point(6) = { widthX/2, -widthY/2, -depthDomain, meshSizeEdges};
Point(7) = {-widthX/2,  widthY/2, -depthDomain, meshSizeEdges};
Point(8) = { widthX/2,  widthY/2, -depthDomain, meshSizeEdges};

// top surface
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// bottom face
Line(5) = {5, 6};
Line(6) = {6, 8};
Line(7) = {8, 7};
Line(8) = {7, 5};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
// vertical edges:
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {4, 8};
Line(12) = {3, 7};

// Face 1
Curve Loop(3) = {6, -11, -2, 10};
Plane Surface(3) = {3};
// Face 2
Curve Loop(4) = {7, -12, -3, 11};
Plane Surface(4) = {4};
// Face 3
Curve Loop(5) = {8, -9, -4, 12};
Plane Surface(5) = {5};
// Face 4
Curve Loop(6) = {5, -10, -1, 9};
Plane Surface(6) = {6};
// Volume:
Surface Loop(2) = {6, 3, 4, 5, 1, 2};
Volume(1) = {2, 3};
// ... end core ............................................

// ... bounding box .................................................

Point(21) = {-boundingBoxX/2, -boundingBoxY/2, 0, meshSizeBoundingBox};
Point(22) = { boundingBoxX/2, -boundingBoxY/2, 0, meshSizeBoundingBox};
Point(23) = {-boundingBoxX/2,  boundingBoxY/2, 0, meshSizeBoundingBox};
Point(24) = { boundingBoxX/2,  boundingBoxY/2, 0, meshSizeBoundingBox};

Point(25) = {-boundingBoxX/2, -boundingBoxY/2, -boundingBoxZ, meshSizeBoundingBox};
Point(26) = { boundingBoxX/2, -boundingBoxY/2, -boundingBoxZ, meshSizeBoundingBox};
Point(27) = {-boundingBoxX/2,  boundingBoxY/2, -boundingBoxZ, meshSizeBoundingBox};
Point(28) = { boundingBoxX/2,  boundingBoxY/2, -boundingBoxZ, meshSizeBoundingBox};

Line(26) = {25, 26};
Line(27) = {26, 28};
Line(28) = {28, 27};
Line(29) = {27, 25};
Line(30) = {21, 22};
Line(32) = {24, 22};
Line(33) = {24, 23};
Line(34) = {23, 21};
Line(35) = {25, 21};
Line(36) = {26, 22};
Line(37) = {27, 23};
Line(38) = {28, 24};

Curve Loop(13) = {26, 36, -30, -35};
Plane Surface(13) = {13};
Curve Loop(14) = {27, 38, 32, -36};
Plane Surface(14) = {14};
Curve Loop(15) = {27, 28, 29, 26};
Plane Surface(15) = {-15};
Curve Loop(16) = {35, -34, -37, 29};
Plane Surface(16) = {16};
Curve Loop(17) = {30, -32, 33, 34};
Plane Surface(17) = {-17, -1};
Curve Loop(18) = {33, -37, -28, 38};
Plane Surface(18) = {-18};
Surface Loop(4) = {4, 2, 6, 3, 5, 13, 15, 14, 18, 17, 16};
Volume(3) = {4};

// Electrodes:

//Physical Surface("top", 13) = {1};
//Physical Surface("faces", 14) = {5, 4, 3, 6, 2};
For k In {0:numElectrodes-1}
    kk = newp;
    Point(kk) = { -lengthLine/2 + k * distanceElectrodes ,  0, 0, meshSizeElectrodes};
    Point{kk} In Surface {1};
EndFor
// ... set tagging ....
Physical Volume("anomaly") = {2};
Physical Volume("domain") = {1, 3};
Physical Surface("faces") = {16, 15, 18, 14};
Physical Surface("surface") = {17, 1};

