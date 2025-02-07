Mesh.MshFileVersion = 2.2;

// ERT setup
ElectrodeSpacing = 5;
NumElectrodes = 40;
LineLength = ElectrodeSpacing * (NumElectrodes-1) ;
// Geometry:
SeamThickness=13;
ExtractionWidth = 200;
ExtractionLength = 200;
RemainderLength = 250;
GoafLength = 100;
ExtractionLength = 200;
RoadWidth = 8;
RoadHeight = 5;
WingRoadLength = 50;
Road2Offset = 90;
ExtraLength = ExtractionLength/2 + ExtractionLength + GoafLength - (ExtractionLength/2 + Road2Offset+RoadWidth) ;
///

//RoadLengthWest = 200;
//RoadLengthEast = RemainderLength +200;
//LineOffset = (RemainderLength-LineLength)/2;
//LineHeight = RoadHeight/2;

//CoreThickness = LineLength  * 0.6;
//CoreWidth = ExtractionWidth + RoadWidth + LineLength * 0.3;
//PaddingWidth = LineLength * 0.8;

// Mesh sizes:
// region of most interest = coal
meshSizeCenter = ElectrodeSpacing/3;
meshSizeCoreClose = ElectrodeSpacing/3;
meshSizeCore = ElectrodeSpacing;

//meshSizeOuterCore = 2 * ElectrodeSpacing;
//meshSizePadding=(RoadLengthWest+RoadWidth+2*PaddingWidth+RoadLengthEast)/20;
//meshSizeElectrodes = ElectrodeSpacing/10;

// Base
Point(1) = {0, -ExtractionWidth/2, 0, meshSizeCenter };
Point(2) = {0, ExtractionWidth/2, 0, meshSizeCenter };
Point(3) = {RemainderLength, -ExtractionWidth/2, 0, meshSizeCenter };
Point(4) = {RemainderLength,  ExtractionWidth/2, 0, meshSizeCenter };
Point(5) = {RemainderLength+ExtractionLength, -ExtractionWidth/2, 0, meshSizeCore };
Point(6) = {RemainderLength+ExtractionLength,  ExtractionWidth/2, 0, meshSizeCore };
Point(7) = {0, -ExtractionWidth/2-RoadWidth, 0, meshSizeCoreClose };
Point(8) = {0, ExtractionWidth/2+RoadWidth, 0, meshSizeCoreClose };
Point(9) = {RemainderLength+ExtractionLength, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(10) = {RemainderLength+ExtractionLength, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
Point(11) = {RemainderLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth, 0, meshSizeCore };
Point(12) = {RemainderLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth, 0, meshSizeCore };
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

Point(31) = {RemainderLength+ExtractionLength+GoafLength, -ExtractionWidth/2-RoadWidth-WingRoadLength, 0, meshSizeCore };
Point(32) = {RemainderLength+ExtractionLength+GoafLength, ExtractionWidth/2+RoadWidth+WingRoadLength, 0, meshSizeCore };