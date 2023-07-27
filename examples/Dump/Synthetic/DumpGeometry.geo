

DumpCrone = 300; // length along the top of the dump
DumpAngleTop=25;
DumpAngleFoot=25;
ResetDumpFromTip=50;
ResetDumpFromFoot=100;

LengthFoot=100;
LengthTop=100;
DepthBase=50;

DumpHeightAboveTopLevel =20;
DumpHeightAboveDumpFoot =100;
CoreWidth=200;

PaddingX=300;
PaddingY=300;
PaddingZ=200;

MeshSizeCore=8;
MeshSizeBBox=50;

NumElectrodes=32;
SurveyLineOffsetFromTopEdge=30;
SurveyLineOffsetFromFootEdge=30;

ElectrodeSpacing=(DumpCrone-SurveyLineOffsetFromTopEdge-SurveyLineOffsetFromFootEdge)/(NumElectrodes-1);
MeshSizeElectrodes= ElectrodeSpacing/10;
If ( MeshSizeElectrodes > MeshSizeCore )
	MeshSizeElectrodes= MeshSizeCore;
EndIf
Printf("Electrode Spacing %g",  ElectrodeSpacing);
Core1=-CoreWidth/2;
Core2=CoreWidth/2;

TopDumpEdgeFromFoot = DumpHeightAboveDumpFoot/Tan(DumpAngleFoot/180*Pi);
TopDumpEdgeFromTop= DumpHeightAboveTopLevel/Tan(DumpAngleTop/180*Pi);
DumpBase= TopDumpEdgeFromFoot + DumpCrone + TopDumpEdgeFromTop;

// FRONT:
Point(1)={0,Core1, 0, MeshSizeCore};
Point(2)={DumpCrone, Core1, 0, MeshSizeCore};
Point(3)={DumpCrone + TopDumpEdgeFromFoot, Core1, -DumpHeightAboveDumpFoot,  MeshSizeCore};
Point(4)={DumpCrone + TopDumpEdgeFromFoot-ResetDumpFromFoot, Core1, -DumpHeightAboveDumpFoot, MeshSizeCore};
Point(5)={-TopDumpEdgeFromTop+ ResetDumpFromTip, Core1, -DumpHeightAboveTopLevel, MeshSizeCore};
Point(6)={-TopDumpEdgeFromTop, Core1 , -DumpHeightAboveTopLevel, MeshSizeCore};

Point(7)={-TopDumpEdgeFromTop-LengthTop, Core1, -DumpHeightAboveTopLevel, MeshSizeCore};
Point(8)={-TopDumpEdgeFromTop-LengthTop, Core1, -DumpHeightAboveDumpFoot-DepthBase, MeshSizeCore};
Point(9)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot, Core1, -DumpHeightAboveDumpFoot-DepthBase, MeshSizeCore};
Point(10)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot, Core1,-DumpHeightAboveDumpFoot, MeshSizeCore};

Line(1) = {8, 9};
Line(2) = {9, 10};
Line(3) = {10, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {6, 5};
Line(10) = {5, 4};
Line(11) = {4, 3};
Curve Loop(1) = {10, 11, 4, 5, 6, 9};
Plane Surface(1) = {1};
Curve Loop(2) = {1, 2, 3, -11, -10, -9, 7, 8};
Plane Surface(2) = {2};



// BACK:
B0=15;
Point(B0+1)={0, Core2, 0, MeshSizeCore};
Point(B0+2)={DumpCrone, Core2,0, MeshSizeCore};
Point(B0+3)={DumpCrone + TopDumpEdgeFromFoot, Core2,-DumpHeightAboveDumpFoot, MeshSizeCore};
Point(B0+4)={DumpCrone + TopDumpEdgeFromFoot-ResetDumpFromFoot, Core2,-DumpHeightAboveDumpFoot, MeshSizeCore};
Point(B0+5)={-TopDumpEdgeFromTop+ ResetDumpFromTip, Core2, -DumpHeightAboveTopLevel, MeshSizeCore};
Point(B0+6)={-TopDumpEdgeFromTop, Core2, -DumpHeightAboveTopLevel, MeshSizeCore};

Point(B0+7)={-TopDumpEdgeFromTop-LengthTop, Core2, -DumpHeightAboveTopLevel, MeshSizeCore};
Point(B0+8)={-TopDumpEdgeFromTop-LengthTop, Core2, -DumpHeightAboveDumpFoot-DepthBase, MeshSizeCore};
Point(B0+9)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot, Core2,-DumpHeightAboveDumpFoot-DepthBase, MeshSizeCore};
Point(B0+10)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot, Core2,-DumpHeightAboveDumpFoot, MeshSizeCore};

Line(B0+1) = {B0+8, B0+9};
Line(B0+2) = {B0+9, B0+10};
Line(B0+3) = {B0+10, B0+3};
Line(B0+4) = {B0+3, B0+2};
Line(B0+5) = {B0+2, B0+1};
Line(B0+6) = {B0+1, B0+6};
Line(B0+7) = {B0+6, B0+7};
Line(B0+8) = {B0+7, B0+8};
Line(B0+9) = {B0+6, B0+5};
Line(B0+10) = {B0+5, B0+4};
Line(B0+11) = {B0+4, B0+3};
Curve Loop(B0+1) = {B0+10, B0+11, B0+4, B0+5, B0+6, B0+9};
Plane Surface(B0+1) = {B0+1};
Curve Loop(B0+2) = {B0+1, B0+2, B0+3, -B0-11, -B0-10, -B0-9, B0+7, B0+8};
Plane Surface(B0+2) = {B0+2};

// Core body:
Line(27) = {23, 8};
Line(28) = {22, 7};
Line(29) = {21, 6};
Line(30) = {16, 1};
Line(31) = {20, 5};
Line(32) = {17, 2};
Line(33) = {19, 4};
Line(34) = {18, 3};
Line(35) = {25, 10};
Line(36) = {24, 9};
Curve Loop(18) = {22, 28, -7, -29};
Plane Surface(18) = {18};
Curve Loop(19) = {29, -6, -30, 21};
Plane Surface(19) = {19};
Curve Loop(20) = {29, 9, -31, -24};
Plane Surface(20) = {20};
Curve Loop(21) = {25, 33, -10, -31};
Plane Surface(21) = {21};
Curve Loop(22) = {30, -5, -32, 20};
Plane Surface(22) = {22};
Curve Loop(23) = {32, -4, -34, 19};
Plane Surface(23) = {23};
Curve Loop(24) = {26, 34, -11, -33};
Plane Surface(24) = {24};
Curve Loop(25) = {18, 34, -3, -35};
Plane Surface(25) = {25};
Curve Loop(26) = {1, -36, -16, 27};
Plane Surface(26) = {26};
Curve Loop(27) = {27, -8, -28, 23};
Plane Surface(27) = {27};
Curve Loop(28) = {35, -2, -36, 17};
Plane Surface(28) = {28};

Surface Loop(3) = {19, 1, 23, 22, 16, 21, 24, 20};
Volume(3) = {3};
Surface Loop(4) = {27, 26, 2, 28, 25, 17, 18, 20, 21, 24};
Volume(4) = {4};
Physical Volume("Base", 37) = {4};
Physical Volume("Dump", 38) = {3};

// Bounding box
B1=40;
Point(B1+1)={0, Core1-PaddingY, 0, MeshSizeBBox};
Point(B1+2)={DumpCrone, Core1-PaddingY, 0, MeshSizeBBox};
Point(B1+3)={DumpCrone + TopDumpEdgeFromFoot, Core1-PaddingY, -DumpHeightAboveDumpFoot, MeshSizeBBox};
Point(B1+7)={-TopDumpEdgeFromTop-LengthTop-PaddingX, Core1-PaddingY, -DumpHeightAboveTopLevel, MeshSizeBBox};
Point(B1+6)={-TopDumpEdgeFromTop, Core1-PaddingY, -DumpHeightAboveTopLevel, MeshSizeBBox};
Point(B1+8)={-TopDumpEdgeFromTop-LengthTop-PaddingX, Core1-PaddingY, -DumpHeightAboveDumpFoot-DepthBase-PaddingZ, MeshSizeBBox};
Point(B1+9)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot+PaddingX, Core1-PaddingY,-DumpHeightAboveDumpFoot-DepthBase-PaddingZ, MeshSizeBBox};
Point(B1+10)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot+PaddingX, Core1-PaddingY,-DumpHeightAboveDumpFoot, MeshSizeBBox};

B2=60;
Point(B2+1)={0, Core2+PaddingY, 0, MeshSizeBBox};
Point(B2+2)={DumpCrone, Core2+PaddingY, 0, MeshSizeBBox};
Point(B2+3)={DumpCrone + TopDumpEdgeFromFoot, Core2+PaddingY, -DumpHeightAboveDumpFoot, MeshSizeBBox};
Point(B2+7)={-TopDumpEdgeFromTop-LengthTop-PaddingX, Core2+PaddingY, -DumpHeightAboveTopLevel, MeshSizeBBox};
Point(B2+6)={-TopDumpEdgeFromTop, Core2+PaddingY, -DumpHeightAboveTopLevel, MeshSizeBBox};
Point(B2+8)={-TopDumpEdgeFromTop-LengthTop-PaddingX, Core2+PaddingY, -DumpHeightAboveDumpFoot-DepthBase-PaddingZ, MeshSizeBBox};
Point(B2+9)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot+PaddingX, Core2+PaddingY,-DumpHeightAboveDumpFoot-DepthBase-PaddingZ, MeshSizeBBox};
Point(B2+10)={DumpCrone + TopDumpEdgeFromFoot+ LengthFoot+PaddingX, Core2+PaddingY,-DumpHeightAboveDumpFoot, MeshSizeBBox};

Line(37) = {69, 68};
Line(38) = {68, 67};
Line(39) = {67, 66};
Line(40) = {66, 61};
Line(41) = {61, 62};
Line(42) = {62, 63};
Line(43) = {63, 70};
Line(44) = {70, 69};
Line(45) = {47, 48};
Line(46) = {48, 49};
Line(47) = {49, 50};
Line(48) = {50, 43};
Line(49) = {43, 42};
Line(50) = {42, 41};
Line(51) = {41, 46};
Line(52) = {46, 47};
Line(53) = {69, 49};
Line(54) = {70, 50};
Line(55) = {63, 18};
Line(56) = {3, 43};
Line(57) = {62, 17};
Line(58) = {2, 42};
Line(59) = {61, 16};
Line(60) = {1, 41};
Line(61) = {66, 21};
Line(62) = {6, 46};
Line(63) = {67, 47};
Line(64) = {68, 48};
Curve Loop(29) = {53, -46, -64, -37};
Plane Surface(29) = {29};
Curve Loop(30) = {64, -45, -63, -38};
Plane Surface(30) = {30};
Curve Loop(31) = {38, 39, 40, 41, 42, 43, 44, 37};
Plane Surface(31) = {31};
Curve Loop(32) = {61, 22, 28, -7, 62, 52, -63, 39};
Plane Surface(32) = {32};
Curve Loop(33) = {61, -21, -59, -40};
Plane Surface(33) = {33};
Curve Loop(34) = {6, 62, -51, -60};
Plane Surface(34) = {34};
Curve Loop(35) = {5, 60, -50, -58};
Plane Surface(35) = {35};
Curve Loop(36) = {57, 20, -59, 41};
Plane Surface(36) = {36};
Curve Loop(37) = {57, -19, -55, -42};
Plane Surface(37) = {37};
Curve Loop(38) = {58, -49, -56, 4};
Plane Surface(38) = {38};
Curve Loop(39) = {48, -56, -3, -35, 18, -55, 43, 54};
Plane Surface(39) = {39};
Curve Loop(40) = {48, 49, 50, 51, 52, 45, 46, 47};
Plane Surface(40) = {40};
Curve Loop(41) = {54, -47, -53, -44};
Plane Surface(41) = {41};
Surface Loop(5) = {41, 39, 40, 38, 35, 34, 32, 33, 36, 37, 31, 30, 29, 28, 2, 26, 17, 27, 1, 16};
Volume(5) = {5};
Physical Surface("OuterFaces", 65) = {29, 31, 32, 40, 41, 30};
Physical Volume("BBox") = {5};

B3=80;
For num In {1:NumElectrodes}
	Point(B3+num)={SurveyLineOffsetFromTopEdge + (num-1) * ElectrodeSpacing, 0, 0, MeshSizeElectrodes};
	Point{B3+num} In Surface{22};
	Physical Point(1000+num) = {B3+num};
EndFor
