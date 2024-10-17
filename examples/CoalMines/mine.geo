Mesh.MshFileVersion = 2.2;

// ERT setup
ElectrodeSpacing = 5;
NumElectrodes = 40;
LineLength = ElectrodeSpacing * (NumElectrodes-1) ;
// Geometry:
ExtractionWidth = 200;
RemainderLength = 220;
RoadWidth = 8;
RoadHeight = 5;
RoadLengthWest = 200;
RoadLengthEast = RemainderLength +200;
LineOffset = 10;
LineHeight = 2;

CoreThickness = LineLength  * 0.8;
CoreWidth = ExtractionWidth + RoadWidth + LineLength * 0.6;
PaddingWidth = LineLength * 0.8;

// Mesh sizes:
//meshSizeCoreCenter = ElectrodeSpacing/3;
//meshSizeCore = ElectrodeSpacing/2;
meshSizeCoreCenter = ElectrodeSpacing;
meshSizeCore = ElectrodeSpacing;
meshSizeOuterCore =  ElectrodeSpacing;
meshSizePadding=(RoadLengthWest+RoadWidth+2*PaddingWidth+RoadLengthEast)/20;
//meshSizeElectrodes = ElectrodeSpacing/10;
//meshSizeElectrodeFace = ElectrodeSpacing/4;

meshSizeElectrodes = ElectrodeSpacing/3;
meshSizeElectrodeFace = ElectrodeSpacing;

// Base
Point(1) = {0, -ExtractionWidth/2, 0, meshSizeElectrodeFace };
Point(2) = {0, ExtractionWidth/2, 0, meshSizeElectrodeFace };
Point(3) = {RemainderLength, -ExtractionWidth/2, 0, meshSizeElectrodeFace };
Point(4) = {RemainderLength,  ExtractionWidth/2, 0, meshSizeElectrodeFace };

Point(5) = {-RoadWidth, -ExtractionWidth/2, 0, meshSizeCore };
Point(6) = {-RoadWidth, ExtractionWidth/2, 0, meshSizeCore };
Point(7) = {-RoadWidth-RoadLengthWest, -ExtractionWidth/2, 0, meshSizeCore };
Point(8) = {-RoadWidth-RoadLengthWest, ExtractionWidth/2, 0, meshSizeCore };

Point(9) = {-RoadWidth-RoadLengthWest, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(10) = {-RoadWidth-RoadLengthWest, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
Point(11) = {RoadLengthEast, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(12) = {RoadLengthEast, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };//+
Line(1) = {1, 3};
Line(2) = {3, 4};
Line(3) = {4, 2};
Line(4) = {2, 1};
Line(5) = {6, 8};
Line(7) = {7, 5};
Line(8) = {5, 6};
Line(9) = {10, 12};
Line(10) = {12, 11};
Line(11) = {11, 9};
Line(12) = {9, 7};
Line(13) = {8, 10};
Curve Loop(37) = {4, 1, 2, 3};
Curve Loop(3) = {9, 10, 11, 12, 7, 8, 5, 13};
Plane Surface(3) = {3, 37};
// Top
Point(21) = {0, -ExtractionWidth/2, RoadHeight, meshSizeElectrodeFace };
Point(22) = {0, ExtractionWidth/2, RoadHeight, meshSizeElectrodeFace };
Point(23) = {RemainderLength, -ExtractionWidth/2, RoadHeight, meshSizeElectrodeFace };
Point(24) = {RemainderLength,  ExtractionWidth/2, RoadHeight, meshSizeElectrodeFace };

Point(25) = {-RoadWidth, -ExtractionWidth/2, RoadHeight, meshSizeCore };
Point(26) = {-RoadWidth, ExtractionWidth/2, RoadHeight, meshSizeCore };
Point(27) = {-RoadWidth-RoadLengthWest, -ExtractionWidth/2, RoadHeight, meshSizeCore };
Point(28) = {-RoadWidth-RoadLengthWest, ExtractionWidth/2, RoadHeight, meshSizeCore };

Point(29) = {-RoadWidth-RoadLengthWest, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCore };
Point(30) = {-RoadWidth-RoadLengthWest, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCore };
Point(31) = {RoadLengthEast, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCore };
Point(32) = {RoadLengthEast, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCore };//+
Line(21) = {21, 23};
Line(22) = {23, 24};
Line(23) = {24, 22};
Line(24) = {22, 21};
Line(25) = {26, 28};
Line(27) = {27, 25};
Line(28) = {25, 26};
Line(29) = {30, 32};
Line(30) = {32, 31};
Line(31) = {31, 29};
Line(32) = {29, 27};
Line(33) = {28, 30};
Curve Loop(23) = {29, 30, 31, 32, 27, 28, 25, 33};
Curve Loop(36) = {22, 23, 24, 21};
Plane Surface(23) = {23,36};


// vertical lines:

Line(34) = {10, 30};
Line(35) = {9, 29};
Line(36) = {12, 32};
Line(37) = {11, 31};
Line(38) = {8, 28};
Line(39) = {7, 27};
Line(40) = {5, 25};
Line(41) = {6, 26};
Line(42) = {2, 22};
Line(43) = {4, 24};
Line(44) = {3, 23};
Line(45) = {1, 21};


Curve Loop(25) = {9, 36, -29, -34};
Plane Surface(25) = {25};
Curve Loop(26) = {10, 37, -30, -36};
Plane Surface(26) = {26};
Curve Loop(27) = {35, -31, -37, 11};
Plane Surface(27) = {27};
Curve Loop(28) = {5, 38, -25, -41};
Plane Surface(28) = {28};
Curve Loop(29) = {8, 41, -28, -40};
Plane Surface(29) = {29};
Curve Loop(30) = {7, 40, -27, -39};
Plane Surface(30) = {30};
Curve Loop(31) = {12, 39, -32, -35};
Plane Surface(31) = {31};
Curve Loop(32) = {13, 34, -33, -38};
Plane Surface(32) = {32};

Curve Loop(33) = {4, 45, -24, -42};
Plane Surface(33) = {33};
Curve Loop(34) = {21, -44, -1, 45};
Plane Surface(34) = {34};
Curve Loop(35) = {2, 43, -22, -44};
Plane Surface(35) = {35};
Curve Loop(38) = {42, -23, -43, 3};
Plane Surface(36) = {38};

// now this is put into a core box
Point(50) = {-RoadLengthWest-RoadWidth, -CoreWidth/2, -CoreThickness/2+RoadHeight/2,meshSizeOuterCore};
Point(51) = {-RoadLengthWest-RoadWidth, CoreWidth/2,  -CoreThickness/2+RoadHeight/2, meshSizeOuterCore};
Point(52) = {RoadLengthEast,  -CoreWidth/2, -CoreThickness/2+RoadHeight/2, meshSizeOuterCore};
Point(53) = {RoadLengthEast,  CoreWidth/2,  -CoreThickness/2+RoadHeight/2, meshSizeOuterCore};
Point(54) = {-RoadLengthWest-RoadWidth, -CoreWidth/2, CoreThickness/2+RoadHeight/2,meshSizeOuterCore};
Point(55) = {-RoadLengthWest-RoadWidth, CoreWidth/2,  CoreThickness/2+RoadHeight/2, meshSizeOuterCore};
Point(56) = {RoadLengthEast,  -CoreWidth/2, CoreThickness/2+RoadHeight/2, meshSizeOuterCore};
Point(57) = {RoadLengthEast,  CoreWidth/2,  CoreThickness/2+RoadHeight/2, meshSizeOuterCore};


// and we add some coarse padding - just to make sure!
Line(46) = {54, 56};
Line(47) = {56, 57};
Line(48) = {57, 55};
Line(49) = {55, 54};
Line(50) = {50, 52};
Line(51) = {52, 53};
Line(52) = {53, 51};
Line(53) = {51, 50};
Line(54) = {50, 54};
Line(55) = {51, 55};
Line(56) = {53, 57};
Line(57) = {52, 56};
Curve Loop(39) = {46, -57, -50, 54};
Plane Surface(37) = {39};
Curve Loop(40) = {52, 55, -48, -56};
Plane Surface(38) = {40};
Curve Loop(41) = {51, 56, -47, -57};
Plane Surface(39) = {41, 26};
Curve Loop(42) = {48, 49, 46, 47};
Plane Surface(40) = {42};
Curve Loop(43) = {50, 51, 52, 53};
Plane Surface(41) = {43};
Curve Loop(44) = {53, 54, -49, -55};
Plane Surface(42) = {44, 31, 32};
Surface Loop(1) = {42, 41, 37, 40, 38, 39, 25, 3, 27, 23, 30, 29, 28, 34, 35, 36, 33};
Volume(1) = {1};

Point(60) = {-RoadLengthWest-RoadWidth-PaddingWidth, -CoreWidth/2-PaddingWidth, -CoreThickness/2+RoadHeight/2-PaddingWidth,meshSizePadding};
Point(61) = {-RoadLengthWest-RoadWidth-PaddingWidth, CoreWidth/2+PaddingWidth,  -CoreThickness/2+RoadHeight/2-PaddingWidth, meshSizePadding};
Point(62) = {RoadLengthEast+PaddingWidth,  -CoreWidth/2-PaddingWidth, -CoreThickness/2+RoadHeight/2-PaddingWidth, meshSizePadding};
Point(63) = {RoadLengthEast+PaddingWidth,  CoreWidth/2+PaddingWidth,  -CoreThickness/2+RoadHeight/2-PaddingWidth, meshSizePadding};
Point(64) = {-RoadLengthWest-RoadWidth-PaddingWidth, -CoreWidth/2-PaddingWidth, CoreThickness/2+RoadHeight/2+PaddingWidth,meshSizePadding};
Point(65) = {-RoadLengthWest-RoadWidth-PaddingWidth, CoreWidth/2+PaddingWidth,  CoreThickness/2+RoadHeight/2+PaddingWidth, meshSizePadding};
Point(66) = {RoadLengthEast+PaddingWidth,  -CoreWidth/2-PaddingWidth, CoreThickness/2+RoadHeight/2+PaddingWidth, meshSizePadding};
Point(67) = {RoadLengthEast+PaddingWidth,  CoreWidth/2+PaddingWidth,  CoreThickness/2+RoadHeight/2+PaddingWidth, meshSizePadding};
//+
Line(58) = {64, 60};
Line(59) = {60, 61};
Line(60) = {61, 65};
Line(61) = {65, 64};
Line(62) = {64, 66};
Line(63) = {66, 67};
Line(64) = {67, 65};
Line(65) = {61, 63};
Line(66) = {63, 67};
Line(67) = {60, 62};
Line(68) = {62, 66};
Line(69) = {63, 62};
Curve Loop(45) = {60, -64, -66, -65};
Plane Surface(43) = {45};
Curve Loop(46) = {61, 58, 59, 60};
Plane Surface(44) = {-46};
Curve Loop(47) = {62, -68, -67, -58};
Plane Surface(45) = {-47};
Curve Loop(48) = {69, 68, 63, -66};
Plane Surface(46) = {-48};
Curve Loop(49) = {69, -67, 59, 65};
Plane Surface(47) = {49};
Curve Loop(50) = {61, 62, 63, 64};
Plane Surface(48) = {50};
Surface Loop(2) = {44, 43, 48, 45, 47, 46};
Surface Loop(3) = {37, 40, 38, 41, 39, 42, 26, 31, 32};
Volume(2) = {2, 3};

Physical Volume("CoreRegion", 70) = {1};
Physical Volume("PaddingRegion", 71) = {2};
Physical Surface("Faces", 72) = {44, 45, 47, 46, 48, 43};

For k In {0:NumElectrodes-1}
    kk = newp;
    Point(kk) = { LineOffset + k * ElectrodeSpacing ,   -ExtractionWidth/2, LineHeight, meshSizeElectrodes};
    Point{kk} In Surface {34};
    kk = newp;
    Point(kk) = { LineOffset + k * ElectrodeSpacing ,   ExtractionWidth/2, LineHeight, meshSizeElectrodes};
    Point{kk} In Surface {36};
EndFor

