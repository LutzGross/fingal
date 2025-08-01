// GMSH geometry file
// 
Mesh.MshFileVersion = 2.2;
// ... electrodes in parallel lines .....
numLines = 3;
distanceLines = 4;
numElectrodes = 32;
distanceElectrodes = 3;
//
lengthLine = (numElectrodes-1) * distanceElectrodes;
widthLines = (numLines -1) * distanceLines;
// ... Anomaly ...
anomalyDepth = 2 * distanceElectrodes;
anomalyOffsetX =  -6 * distanceElectrodes;
anomalyLengthX = 16 * distanceElectrodes;
anomalyLengthY = 12 * distanceElectrodes + widthLines;
anomalyLengthZ = 4 * distanceElectrodes;

// ... core region ....
depthDomain  = lengthLine * 0.30;
widthY = anomalyLengthY + 2 * distanceElectrodes;
widthX = lengthLine + 8 * distanceElectrodes;

// ... Padding ......
paddingWidth = 75 * distanceElectrodes ;
boundingBoxX =  widthX + 2 * paddingWidth;
boundingBoxY =  widthY + 2 * paddingWidth;
boundingBoxZ =  depthDomain + paddingWidth;

// ... element sizes ...
meshSizeEdges = distanceElectrodes * 0.5;
meshSizeElectrodes = distanceElectrodes * 0.1;
meshSizeBoundingBox = boundingBoxY/10;

// .... hexahedral anomaly right ........................................................................
Point(109) = { anomalyOffsetX-anomalyLengthX/2, -anomalyLengthY/2, -anomalyDepth, 		meshSizeEdges};
Point(110) = { anomalyOffsetX+ anomalyLengthX/2, -anomalyLengthY/2, -anomalyDepth, meshSizeEdges};
Point(111) = {anomalyOffsetX-anomalyLengthX/2,  anomalyLengthY/2, -anomalyDepth, 		meshSizeEdges};
Point(112) = { anomalyOffsetX+ anomalyLengthX/2,  anomalyLengthY/2, -anomalyDepth, meshSizeEdges};

Point(113) = {anomalyOffsetX-anomalyLengthX/2, -anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 		meshSizeEdges};
Point(114) = { anomalyOffsetX+ anomalyLengthX/2, -anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 	meshSizeEdges};
Point(115) = {anomalyOffsetX-anomalyLengthX/2,  anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 		meshSizeEdges};
Point(116) = { anomalyOffsetX+ anomalyLengthX/2,  anomalyLengthY/2, -anomalyDepth-anomalyLengthZ, 	meshSizeEdges};

Line(113) = {109, 110};
Line(114) = {110, 112};
Line(115) = {112, 111};
Line(116) = {111, 109};
Line(117) = {113, 114};
Line(118) = {113, 114};
Line(119) = {114, 116};
Line(120) = {116, 115};
Line(121) = {115, 113};
Line(122) = {113, 109};
Line(123) = {114, 110};
Line(124) = {116, 112};
Line(125) = {115, 111};
Curve Loop(107) = {117, 123, -113, -122};
Plane Surface(107) = {107};
Curve Loop(108) = {119, 124, -114, -123};
Plane Surface(108) = {108};
Curve Loop(109) = {120, 125, -115, -124};
Plane Surface(109) = {109};
Curve Loop(110) = {116, 113, 114, 115};
Plane Surface(110) = {110};
Curve Loop(111) = {125, 116, -122, -121};
Plane Surface(111) = {111};
Curve Loop(112) = {117, 119, 120, 121};
Plane Surface(112) = {112};
Surface Loop(103) = {111, 109, 112, 107, 108, 110};
Volume(102) = {103};
// ........ End Anomaly Left................

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
Volume(1) = {2, 103};
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
For l In {0:numLines-1}
    For k In {0:numElectrodes-1}
        kk = newp;
        Point(kk) = { -lengthLine/2 + k * distanceElectrodes ,   -widthLines/2 + l * distanceLines, 0, meshSizeElectrodes};
        Point{kk} In Surface {1};
    EndFor
EndFor
// ... set tagging ....
Physical Volume("anomaly") = {102};
Physical Volume("domain") = {1};
Physical Volume("padding") = {3};
Physical Surface("faces") = {16, 15, 18, 14};
Physical Surface("surface") = {17, 1};

