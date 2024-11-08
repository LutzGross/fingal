Mesh.MshFileVersion = 2.2;
L = 300;
Depth = 200;

L_anomaly = 24;
Depth_anomaly = 5;
Offset_anomaly = 0;

LineOffset = 0;
NumElectrodes =  8;
ElectrodeSpacing = 5;
LineLength = (NumElectrodes-1) * ElectrodeSpacing;

meshSizeEdges = L/10;
meshSizeElectrodes = ElectrodeSpacing/30;
meshSizeEdges2 = L_anomaly/10;

Point(1) = {-L/2, L/2, 0, meshSizeEdges};
Point(2) = {-L/2, -L/2, 0, meshSizeEdges};
Point(3) = {L/2, -L/2, 0, meshSizeEdges};
Point(4) = {L/2, L/2, 0, meshSizeEdges};
Point(11) = {-L/2, L/2, -Depth, meshSizeEdges};
Point(12) = {-L/2, -L/2,-Depth, meshSizeEdges};
Point(13) = {L/2, -L/2, -Depth, meshSizeEdges};
Point(14) = {L/2, L/2, -Depth, meshSizeEdges};

Point(31) = {-L_anomaly/2+Offset_anomaly, L/3, -Depth_anomaly, meshSizeEdges2};
Point(32) = {-L_anomaly/2+Offset_anomaly, -L/3, -Depth_anomaly, meshSizeEdges2};
Point(33) = {L_anomaly/2+Offset_anomaly, -L/3, -Depth_anomaly, meshSizeEdges2};
Point(34) = {L_anomaly/2+Offset_anomaly, L/3, -Depth_anomaly, meshSizeEdges2};
Point(41) = {-L_anomaly/2+Offset_anomaly, L/3, -Depth_anomaly-L_anomaly, meshSizeEdges2};
Point(42) = {-L_anomaly/2+Offset_anomaly, -L/3, -Depth_anomaly-L_anomaly, meshSizeEdges2};
Point(43) = {L_anomaly/2+Offset_anomaly, -L/3, -Depth_anomaly-L_anomaly, meshSizeEdges2};
Point(44) = {L_anomaly/2+Offset_anomaly, L/3, -Depth_anomaly-L_anomaly, meshSizeEdges2};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

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
Surface Loop(1) = {14, 12, 11, 15, 1, 13};

Line(19) = {42, 43};
Line(20) = {43, 33};
Line(21) = {33, 32};
Line(22) = {32, 42};
Line(23) = {34, 33};
Line(24) = {44, 43};
Line(25) = {41, 42};
Line(26) = {32, 31};
Line(27) = {31, 34};
Line(28) = {34, 44};
Line(29) = {31, 41};
Line(30) = {41, 44};
Curve Loop(16) = {25, -22, 26, 29};
Plane Surface(16) = {16};
Curve Loop(17) = {30, -28, -27, 29};
Plane Surface(17) = {17};
Curve Loop(18) = {24, 20, -23, 28};
Plane Surface(18) = {18};
Curve Loop(19) = {19, 20, 21, 22};
Plane Surface(19) = {19};
Curve Loop(20) = {21, 26, 27, 23};
Plane Surface(20) = {20};
Curve Loop(21) = {19, -24, -30, 25};
Plane Surface(21) = {21};
Surface Loop(2) = {16, 21, 19, 18, 20, 17};
Volume(2) = {2};
Volume(1) = {1,2};

For k In {0:NumElectrodes-1}
    kk = newp;
    Point(kk) = { -LineLength/2 + k * ElectrodeSpacing ,   LineOffset, 0, meshSizeElectrodes};
    Point{kk} In Surface {1};
EndFor

Physical Surface("faces", 19) = {13, 11, 12, 14, 15};
Physical Volume("core", 20) = {1};
Physical Volume("anomaly", 21) = {2};
