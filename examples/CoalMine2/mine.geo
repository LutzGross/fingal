Mesh.MshFileVersion = 2.2;

// ERT setup
ElectrodeSpacing = 5;
NumElectrodes = 49;
LineOffset = 5;
LineHeight = 1;
LineLength = ElectrodeSpacing * (NumElectrodes-1) ;
// Geometry:
SeamThickness=13;
ExtractionWidth = 200;
ExtractionLength = 10;
UnMinedLength = 250;
GoafLength = 250;
RoadWidth = 5;
RoadHeight = 5;
WingRoadLength = 50;
Road2Offset = 90;
BaseThickness = 100;
RockMassThickness = 100;
PaddingWidth = 100;
// not used for meshing but read to generate synthetic data:
DamageZoneOffset = 150;
DamageHeight = 80;
DamageBaseDepth = 40;
DamageGeometryExponent = 0.5;
FringeWidth = 5;
ResetDamagedZoneSouth = 3;
ResetDamagedZoneNorth =3;

ExtraLength = ExtractionLength/2 + ExtractionLength + GoafLength - (ExtractionLength/2 + Road2Offset+RoadWidth) ;
CoreLengthX=UnMinedLength+ExtractionLength+GoafLength + Road2Offset+RoadWidth+ExtraLength;
CoreLengthY=ExtractionWidth+2*RoadWidth+2*WingRoadLength;
// Mesh sizes:
// region of most interest = coal

meshSizeCenter = ElectrodeSpacing;
meshSizeCenterNorth = ElectrodeSpacing/5;
meshSizeCenterSouth = ElectrodeSpacing/5;
meshSizeCoreClose = ElectrodeSpacing;
meshSizeCore = 2*ElectrodeSpacing;
meshSizePadding=(CoreLengthX+2*PaddingWidth)/15;
meshSizeElectrodes = ElectrodeSpacing/25;

// Base
Point(1) = {0, -ExtractionWidth/2, 0, meshSizeCenterSouth };
Point(2) = {0, ExtractionWidth/2, 0, meshSizeCenterNorth };
Point(3) = {UnMinedLength, -ExtractionWidth/2, 0, meshSizeCenterSouth };
Point(4) = {UnMinedLength,  ExtractionWidth/2, 0, meshSizeCenterNorth };
Point(7) = {0, -ExtractionWidth/2-RoadWidth, 0, meshSizeCoreClose };
Point(8) = {0, ExtractionWidth/2+RoadWidth, 0, meshSizeCoreClose };
Point(9) = {UnMinedLength+ExtractionLength, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(10) = {UnMinedLength+ExtractionLength, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
Point(11) = {UnMinedLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(12) = {UnMinedLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
Point(13) = {0, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(14) = {0, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };
Point(15) = {-RoadWidth, -ExtractionWidth/2-RoadWidth, 0, meshSizeCoreClose };
Point(16) = {-RoadWidth, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };


Point(17) = {-RoadWidth, -ExtractionWidth/2, 0, meshSizeCoreClose };
Point(18) = {-RoadWidth, ExtractionWidth/2, 0, meshSizeCoreClose };
Point(19) = {-RoadWidth, ExtractionWidth/2+RoadWidth, 0, meshSizeCoreClose };
Point(20) = {-RoadWidth, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };

Point(21) = {-Road2Offset, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(22) = {-Road2Offset, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(23) = {-Road2Offset, -ExtractionWidth/2, 0, meshSizeCore };
Point(24) = {-Road2Offset, ExtractionWidth/2, 0, meshSizeCore };
Point(25) = {-Road2Offset, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
Point(26) = {-Road2Offset, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };

Point(27) = {-Road2Offset-RoadWidth, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(28) = {-Road2Offset-RoadWidth, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };

Point(29) = {-Road2Offset-RoadWidth-ExtraLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(30) = {-Road2Offset-RoadWidth-ExtraLength, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };

Point(31) = {UnMinedLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(32) = {UnMinedLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };

Point(34) = {UnMinedLength, ExtractionWidth/2+RoadWidth, 0, meshSizeCoreClose };
Point(33) = {UnMinedLength, -ExtractionWidth/2-RoadWidth, 0, meshSizeCoreClose };

Line(1) = {30, 28};
Line(2) = {28, 26};
Line(3) = {26, 20};
Line(4) = {20, 14};
Line(5) = {14, 32};
Line(6) = {25, 19};
Line(7) = {8, 34};
Line(40) = {34, 10};
Line(9) = {10, 12};
Line(10) = {24, 18};
Line(12) = {23, 17};
Line(13) = {1, 3};
Line(15) = {22, 15};
Line(16) = {7, 33};
Line(41) = {33, 9};
Line(17) = {9, 11};
Line(18) = {29, 27};
Line(19) = {27, 21};
Line(20) = {21, 16};
Line(21) = {16, 13};
Line(22) = {13, 31};
Line(23) = {30, 29};
Line(24) = {28, 27};
Line(25) = {26, 25};
Line(27) = {24, 23};
Line(28) = {22, 21};
Line(29) = {20, 19};
Line(30) = {18, 17};
Line(31) = {15, 16};
Line(32) = {14, 8};
Line(33) = {2, 1};
Line(34) = {7, 13};
Line(35) = {4, 3};
Line(36) = {10, 9};
Line(37) = {32, 12};
Line(38) = {12, 11};
Line(39) = {11, 31};
Line(42) = {2, 4};
Line(43) = {4, 34};
Line(44) = {33, 3};
Curve Loop(1) = {24, -18, -23, 1};
Plane Surface(1) = {1};
Curve Loop(2) = {24, 19, -28, 15, 31, 21, -34, 16, 44, -13, -33, 42, 43, -7, -32, -4, 29, -6, -25, -2};
Curve Loop(3) = {10, 30, -12, -27};
Plane Surface(2) = {2, 3};
Plane Surface(3) = {3};
Curve Loop(4) = {6, -29, -3, 25};
Plane Surface(4) = {4};
Curve Loop(5) = {5, 37, -9, -40, -7, -32};
Plane Surface(5) = {5};
Curve Loop(6) = {40, 36, -41, 44, -35, 43};
Plane Surface(6) = {6};
Curve Loop(7) = {35, -13, -33, 42};
Plane Surface(7) = {7};
Curve Loop(8) = {16, 41, 17, 39, -22, -34};
Plane Surface(8) = {8};
Curve Loop(9) = {15, 31, -20, -28};
Plane Surface(9) = {9};
Curve Loop(10) = {17, -38, -9, 36};
Plane Surface(10) = {10};

// Roads
Point(101) = {0, -ExtractionWidth/2, RoadHeight, meshSizeCenterSouth };
Point(102) = {0, ExtractionWidth/2, RoadHeight, meshSizeCenterNorth };
Point(103) = {UnMinedLength, -ExtractionWidth/2, RoadHeight, meshSizeCenterSouth };
Point(104) = {UnMinedLength,  ExtractionWidth/2, RoadHeight, meshSizeCenterNorth };
Point(107) = {0, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCoreClose };
Point(108) = {0, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCoreClose };
Point(113) = {0, -ExtractionWidth/2-RoadWidth-WingRoadLength, RoadHeight, meshSizeCore };
Point(114) = {0, ExtractionWidth/2+RoadWidth+WingRoadLength, RoadHeight, meshSizeCore };
Point(115) = {-RoadWidth, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCoreClose };
Point(116) = {-RoadWidth, -ExtractionWidth/2-RoadWidth-WingRoadLength, RoadHeight, meshSizeCore };
Point(117) = {-RoadWidth, -ExtractionWidth/2, RoadHeight, meshSizeCoreClose };
Point(118) = {-RoadWidth, ExtractionWidth/2, RoadHeight, meshSizeCoreClose };
Point(119) = {-RoadWidth, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCoreClose };
Point(120) = {-RoadWidth, ExtractionWidth/2+RoadWidth+WingRoadLength, RoadHeight, meshSizeCore };
Point(121) = {-Road2Offset, -ExtractionWidth/2-RoadWidth-WingRoadLength, RoadHeight, meshSizeCore };
Point(122) = {-Road2Offset, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCore };
Point(123) = {-Road2Offset, -ExtractionWidth/2, RoadHeight, meshSizeCore };
Point(124) = {-Road2Offset, ExtractionWidth/2, RoadHeight, meshSizeCore };
Point(125) = {-Road2Offset, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCore };
Point(126) = {-Road2Offset, ExtractionWidth/2+RoadWidth+WingRoadLength, RoadHeight, meshSizeCore };
Point(127) = {-Road2Offset-RoadWidth, -ExtractionWidth/2-RoadWidth-WingRoadLength, RoadHeight, meshSizeCore };
Point(128) = {-Road2Offset-RoadWidth, ExtractionWidth/2+RoadWidth+WingRoadLength, RoadHeight, meshSizeCore };
Point(134) = {UnMinedLength, ExtractionWidth/2+RoadWidth, RoadHeight, meshSizeCoreClose };
Point(133) = {UnMinedLength, -ExtractionWidth/2-RoadWidth, RoadHeight, meshSizeCoreClose };//+
Line(45) = {128, 126};
Line(46) = {126, 125};
Line(47) = {125, 119};
Line(48) = {119, 120};
Line(49) = {120, 114};
Line(50) = {114, 108};
Line(51) = {108, 134};
Line(52) = {134, 104};
Line(53) = {104, 102};
Line(54) = {102, 101};
Line(55) = {101, 103};
Line(56) = {103, 133};
Line(57) = {133, 107};
Line(58) = {107, 113};
Line(59) = {113, 116};
Line(60) = {116, 115};
Line(61) = {115, 122};
Line(62) = {122, 121};
Line(63) = {121, 127};
Line(64) = {127, 128};
Line(65) = {128, 28};
Line(66) = {126, 26};
Line(67) = {125, 25};
Line(68) = {24, 124};
Line(69) = {23, 123};
Line(70) = {123, 117};
Line(71) = {117, 17};
Line(72) = {117, 118};
Line(73) = {118, 18};
Line(74) = {118, 124};
Line(75) = {124, 123};
Line(76) = {8, 108};
Line(77) = {2, 102};
Line(78) = {134, 34};
Line(79) = {104, 4};
Line(80) = {14, 114};
Line(81) = {20, 120};
Line(82) = {1, 101};
Line(83) = {7, 107};
Line(84) = {13, 113};
Line(85) = {16, 116};
Line(86) = {21, 121};
Line(87) = {27, 127};
Line(88) = {22, 122};
Line(89) = {15, 115};
Line(90) = {33, 133};
Line(91) = {3, 103};
Line(92) = {19, 119};

Curve Loop(11) = {86, 63, -87, 19};
Plane Surface(11) = {11};
Curve Loop(12) = {28, 86, -62, -88};
Plane Surface(12) = {12};
Curve Loop(13) = {15, 89, 61, -88};
Plane Surface(13) = {13};
Curve Loop(14) = {89, -60, -85, -31};
Plane Surface(14) = {14};
Curve Loop(15) = {59, -85, 21, 84};
Plane Surface(15) = {15};
Curve Loop(16) = {34, 84, -58, -83};
Plane Surface(16) = {16};
Curve Loop(17) = {16, 90, 57, -83};
Plane Surface(17) = {17};
//Curve Loop(18) = {44, 91, 56, -90};
//Plane Surface(18) = {18};
Curve Loop(19) = {13, 91, -55, -82};
Plane Surface(19) = {19};
Curve Loop(20) = {33, 82, -54, -77};
Plane Surface(20) = {20};
Curve Loop(21) = {53, -77, 42, -79};
Plane Surface(21) = {21};
//Curve Loop(22) = {79, 43, -78, 52};
//Plane Surface(22) = {22};
Curve Loop(23) = {7, -78, -51, -76};
Plane Surface(23) = {23};
Curve Loop(24) = {32, 76, -50, -80};
Plane Surface(24) = {24};
Curve Loop(25) = {80, -49, -81, 4};
Plane Surface(25) = {25};
Curve Loop(26) = {48, -81, 29, 92};
Plane Surface(26) = {26};
Curve Loop(27) = {6, 92, -47, 67};
Plane Surface(27) = {27};
Curve Loop(28) = {25, -67, -46, 66};
Plane Surface(28) = {28};
Curve Loop(29) = {2, -66, -45, 65};
Plane Surface(29) = {29};
Curve Loop(30) = {12, -71, -70, -69};
Plane Surface(30) = {30};
Curve Loop(31) = {27, 69, -75, -68};
Plane Surface(31) = {31};
Curve Loop(32) = {10, -73, 74, -68};
Plane Surface(32) = {32};
Curve Loop(33) = {30, -71, 72, 73};
Plane Surface(33) = {33};
Curve Loop(34) = {24, 87, 64, 65};
Plane Surface(34) = {34};

Curve Loop(35) = {45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64};
Curve Loop(36) = {75, 70, 72, 74};
Plane Surface(35) = {35, 36};
//Plane Surface(36) = {36};
// Coal seam:
Point(151) = {0, -ExtractionWidth/2, SeamThickness, meshSizeCenter };
Point(152) = {0, ExtractionWidth/2, SeamThickness, meshSizeCenter };
Point(153) = {UnMinedLength, -ExtractionWidth/2-RoadWidth, SeamThickness, meshSizeCenter };
Point(154) = {UnMinedLength,  ExtractionWidth/2+RoadWidth, SeamThickness, meshSizeCenter };
Point(155) = {UnMinedLength, -ExtractionWidth/2, SeamThickness, meshSizeCenter };
Point(156) = {UnMinedLength,  ExtractionWidth/2, SeamThickness, meshSizeCenter };
Point(159) = {UnMinedLength+ExtractionLength, -ExtractionWidth/2-RoadWidth, SeamThickness, meshSizeCore };
Point(160) = {UnMinedLength+ExtractionLength, ExtractionWidth/2+RoadWidth, SeamThickness, meshSizeCore };
Point(161) = {UnMinedLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth, SeamThickness, meshSizeCore };
Point(162) = {UnMinedLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth, SeamThickness, meshSizeCore };//+
Line(93) = {101, 151};
Line(94) = {102, 152};
Line(95) = {152, 151};
Line(96) = {133, 153};
Line(97) = {153, 159};
Line(98) = {159, 9};
Line(99) = {159, 161};
Line(100) = {161, 11};
Line(101) = {161, 162};
Line(102) = {162, 12};
Line(103) = {162, 160};
Line(104) = {160, 10};
Line(105) = {160, 159};
Line(106) = {160, 154};
Line(107) = {151, 155};
Line(108) = {155, 153};
Line(109) = {154, 156};
Line(110) = {156, 152};
Line(111) = {104, 156};
Line(112) = {103, 155};
Line(113) = {156, 155};
Line(114) = {134, 154};
Curve Loop(37) = {95, -93, -54, 94};
Plane Surface(36) = {37};
Curve Loop(38) = {93, 107, -112, -55};
Plane Surface(37) = {38};
Curve Loop(39) = {112, 108, -96, -56};
Plane Surface(38) = {39};
Curve Loop(40) = {97, 98, -41, 90, 96};
Plane Surface(39) = {40};
Curve Loop(41) = {98, 17, -100, -99};
Plane Surface(40) = {41};
Curve Loop(42) = {36, -98, -105, 104};
Plane Surface(41) = {42};
Curve Loop(43) = {38, -100, 101, 102};
Plane Surface(42) = {43};
Curve Loop(44) = {104, 9, -102, 103};
Plane Surface(43) = {44};
Curve Loop(45) = {105, 99, 101, 103};
Plane Surface(44) = {45};
Curve Loop(46) = {40, -104, 106, -114, 78};
Plane Surface(45) = {46};
Curve Loop(47) = {52, 111, -109, -114};
Plane Surface(46) = {47};
Curve Loop(48) = {79, 35, 91, 112, -113, -111};
Plane Surface(47) = {48};
Curve Loop(49) = {113, -107, -95, -110};
Plane Surface(48) = {49};
Curve Loop(50) = {113, 108, 97, -105, 106, 109};
Plane Surface(49) = {50};
Surface Loop(1) = {42, 10, 40, 44, 43, 41};
// ... Goaf
Volume(1) = {1};
Physical Volume("Goaf", 115) = {1};
// Seam
Curve Loop(51) = {53, 94, -110, -111};
Plane Surface(50) = {51};
Surface Loop(2) = {48, 37, 36, 50, 47, 21, 20, 19, 7};
Volume(2) = {2};
Physical Volume("Seam", 116) = {2};
// base
Point(329) = {-Road2Offset-RoadWidth-ExtraLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, -BaseThickness, meshSizeCore };
Point(330) = {-Road2Offset-RoadWidth-ExtraLength, ExtractionWidth/2+RoadWidth+WingRoadLength, -BaseThickness, meshSizeCore };
Point(331) = {UnMinedLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, -BaseThickness, meshSizeCore };
Point(332) = {UnMinedLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth+WingRoadLength, -BaseThickness, meshSizeCore };//+
Line(115) = {329, 331};
Line(116) = {331, 332};
Line(117) = {332, 330};
Line(118) = {330, 329};
Line(119) = {329, 29};
Line(120) = {330, 30};
Line(121) = {332, 32};
Line(122) = {331, 31};
Curve Loop(52) = {115, 122, -22, -21, -20, -19, -18, -119};
Plane Surface(51) = {52};
Curve Loop(53) = {119, -23, -120, 118};
Plane Surface(52) = {53};
Curve Loop(54) = {5, -121, 117, 120, 1, 2, 3, 4};
Plane Surface(53) = {54};
Curve Loop(55) = {116, 121, 37, 38, 39, -122};
Plane Surface(54) = {55};
Curve Loop(56) = {115, 116, 117, 118};
Plane Surface(55) = {56};
Surface Loop(3) = {-1,2,-3,4,-5,-6,-7,-8,-9,10, 52,53,54,-55,51};
Volume(3) = {3};
Physical Volume("Base", 117) = {3};
// MAss
Point(340) = {-Road2Offset-RoadWidth-ExtraLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, RockMassThickness, meshSizeCore };
Point(341) = {-Road2Offset-RoadWidth-ExtraLength, ExtractionWidth/2+RoadWidth+WingRoadLength, RockMassThickness, meshSizeCore };
Point(342) = {UnMinedLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, RockMassThickness, meshSizeCore };
Point(343) = {UnMinedLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth+WingRoadLength, RockMassThickness, meshSizeCore };
Line(123) = {29, 340};
Line(124) = {30, 341};
Line(125) = {341, 343};
Line(126) = {343, 32};
Line(127) = {343, 342};
Line(128) = {342, 31};
Line(129) = {342, 340};
Line(130) = {340, 341};
Curve Loop(57) = {130, -124, 23, 123};
Plane Surface(56) = {57};
Curve Loop(58) = {1, -65, 45, 66, 3, 81, 49, -80, 5, -126, -125, -124};
Plane Surface(57) = {58};
Curve Loop(59) = {37, -102, -101, 100, 39, -128, -127, 126};
Plane Surface(58) = {59};
Curve Loop(60) = {128, -22, 84, 59, -85, -20, 86, 63, -87, -18, 123, -129};
Plane Surface(59) = {60};
Curve Loop(61) = {130, 125, 127, 129};
Plane Surface(60) = {61};
Surface Loop(4) ={-57, 3, 1, 48, 9, 14, 8, -17, 16, -13, 12, -33, 30, 31, -37, -36, -26, 27, 4, 28, -34, 5, 49, 39, -40, -32, 35, -44, 38, 23, -46, -43, 45, -50, 24, -59, -60, 58, 56};
Volume(4) = {4};
Physical Volume("Mass", 118) = {4};

// Padding


Point(441) = {-Road2Offset-RoadWidth-ExtraLength-PaddingWidth, -CoreLengthY/2-PaddingWidth, -BaseThickness-PaddingWidth, meshSizePadding };
Point(442) = { UnMinedLength+ExtractionLength+GoafLength+PaddingWidth, -CoreLengthY/2-PaddingWidth, -BaseThickness-PaddingWidth, meshSizePadding };
Point(443) = {-Road2Offset-RoadWidth-ExtraLength-PaddingWidth, CoreLengthY/2+PaddingWidth, -BaseThickness-PaddingWidth, meshSizePadding };
Point(444) = { UnMinedLength+ExtractionLength+GoafLength+PaddingWidth, CoreLengthY/2+PaddingWidth, -BaseThickness-PaddingWidth, meshSizePadding };
Point(445) = {-Road2Offset-RoadWidth-ExtraLength-PaddingWidth, -CoreLengthY/2-PaddingWidth, RockMassThickness+PaddingWidth, meshSizePadding };
Point(446) = { UnMinedLength+ExtractionLength+GoafLength+PaddingWidth, -CoreLengthY/2-PaddingWidth, RockMassThickness+PaddingWidth, meshSizePadding };
Point(447) = {-Road2Offset-RoadWidth-ExtraLength-PaddingWidth, CoreLengthY/2+PaddingWidth, RockMassThickness+PaddingWidth, meshSizePadding };
Point(448) = { UnMinedLength+ExtractionLength+GoafLength+PaddingWidth, CoreLengthY/2+PaddingWidth, RockMassThickness+PaddingWidth, meshSizePadding };

Line(131) = {447, 448};
Line(132) = {448, 446};
Line(133) = {446, 445};
Line(134) = {445, 447};
Line(135) = {447, 443};
Line(136) = {448, 444};
Line(137) = {446, 442};
Line(138) = {445, 441};
Line(139) = {441, 442};
Line(140) = {441, 443};
Line(141) = {443, 444};
Line(142) = {444, 442};
Curve Loop(62) = {131, 132, 133, 134};
Plane Surface(61) = {-62};
Curve Loop(63) = {132, 137, -142, -136};
Plane Surface(62) = {63};
Curve Loop(64) = {136, -141, -135, 131};
Plane Surface(63) = {64};
Curve Loop(65) = {140, 141, 142, -139};
Plane Surface(64) = {65};
Curve Loop(66) = {135, -140, -138, 134};
Plane Surface(65) = {66};
Curve Loop(67) = {133, 138, 139, -137};
Plane Surface(66) = {67};
Surface Loop(5) = {61, 63, 62, 66, 65, 64};
Surface Loop(6) = {-56,-52, -51, 59, -11, -15, 29, 25, 57, -53, 60, 55, -54, 58, 42};
Volume(6) = {5, 6};
Physical Volume("Padding", 119) = {6};
Physical Surface("Faces", 143) = {65, 61, 63, 62, 64, 66};

For k In {0:NumElectrodes-1}
    kk = newp;
    Point(kk) = { LineOffset + k * ElectrodeSpacing ,   -ExtractionWidth/2, LineHeight, meshSizeElectrodes};
    Point{kk} In Surface {19};
EndFor

For k In {0:NumElectrodes-1}
    kk = newp;
    Point(kk) = { LineOffset + k * ElectrodeSpacing ,   ExtractionWidth/2, LineHeight, meshSizeElectrodes};
    Point{kk} In Surface {21};
EndFor
