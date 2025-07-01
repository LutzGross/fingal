Mesh.MshFileVersion = 2.2;

%HEADER%

//outer box points
Point(1) = {XminBB, YminBB, ZCore-DepthBoundingBox, meshSizeBB};
Point(2) = { XmaxBB, YminBB, ZCore-DepthBoundingBox, meshSizeBB};
Point(3) = {XminBB,  YmaxBB, ZCore-DepthBoundingBox, meshSizeBB};
Point(4) = { XmaxBB,  YmaxBB, ZCore-DepthBoundingBox, meshSizeBB};

Point(5) = {XminBB, YminBB, ZCore, meshSizeBB};
Point(6) = { XmaxBB, YminBB, ZCore, meshSizeBB};
Point(7) = {XminBB,  YmaxBB, ZCore, meshSizeBB};
Point(8) = { XmaxBB,  YmaxBB, ZCore, meshSizeBB};

//inner box points
Point(101) = {XminCore, YminCore, ZCore-CoreThickness, meshSizeCore};
Point(102) = { XmaxCore, YminCore, ZCore-CoreThickness, meshSizeCore};
Point(103) = {XminCore,  YmaxCore, ZCore-CoreThickness, meshSizeCore};
Point(104) = { XmaxCore,  YmaxCore, ZCore-CoreThickness, meshSizeCore};

Point(105) = {XminCore, YminCore, ZCore, meshSizeCore};
Point(106) = { XmaxCore, YminCore, ZCore, meshSizeCore};
Point(107) = {XminCore,  YmaxCore, ZCore, meshSizeCore};
Point(108) = { XmaxCore,  YmaxCore, ZCore, meshSizeCore};

//inner box points
Point(201) = {XminCore, YminCore, ZCore-DepthInterface, meshSizeCore};
Point(202) = { XmaxCore, YminCore, ZCore-DepthInterface, meshSizeCore};
Point(203) = {XminCore,  YmaxCore, ZCore-DepthInterface, meshSizeCore};
Point(204) = { XmaxCore,  YmaxCore, ZCore-DepthInterface, meshSizeCore};



//inner conduite at interface:

Point(301) = {OffsetXConduit,       OffsetYConduit,ZCore-DepthInterface, meshSizeConduit};
Point(302) = {OffsetXConduit+RadiusConduit,OffsetYConduit,ZCore-DepthInterface, meshSizeConduit};
Point(303) = {OffsetXConduit,       OffsetYConduit+RadiusConduit,ZCore-DepthInterface, meshSizeConduit};
Point(304) = {OffsetXConduit-RadiusConduit,OffsetYConduit,       ZCore-DepthInterface, meshSizeConduit};
Point(305) = {OffsetXConduit,       OffsetYConduit-RadiusConduit,ZCore-DepthInterface, meshSizeConduit};

Circle(301) = {302,301,303};
Circle(302) = {303,301,304};
Circle(303) = {304,301,305};
Circle(304) = {305,301,302};

Line Loop(305) = {301,302,303,304};
Plane Surface(306) = {-305};

Extrude {0,0,DepthInterface-CoreThickness} {
    Surface{306};
    }
Extrude {0, 0, CoreThickness-DepthBoundingBox} {
    Surface{328}; 
    }
//Characteristic Length { 326, 331, 319, 321, 326 } = RadiusConduit;
//Characteristic Length { 321, 326 } = meshSizeCore;
//Characteristic Length { 331, 319  } = RadiusConduit;

//make outside box lines
//bottom
Line(1) = {1,2} ;
Line(2) = {3,4} ;
Line(3) = {1,3} ;
Line(4) = {2,4} ;
//top
Line(5) = {5,6} ;
Line(6) = {7,8} ;
Line(7) = {5,7} ;
Line(8) = {6,8} ;
//sides
Line(9) = {1,5} ;
Line(10) = {2,6} ;
Line(11) = {3,7} ;
Line(12) = {4,8} ;

// outward normals!!! bottom top front back left right
Line Loop(13) = {3,2,-4,-1} ;   
Line Loop(14) = {5,8,-6,-7} ;   
Line Loop(15) = {1,10,-5,-9} ;  
Line Loop(16) = {11,6,-12,-2} ; 
Line Loop(17) = {9,7,-11,-3} ;  
Line Loop(18) = {4,12,-8,-10} ; 



Line(345) = {108, 106};
Line(346) = {204, 202};
Line(347) = {104, 102};
Line(348) = {106, 105};
Line(349) = {105, 107};
Line(350) = {107, 108};
Line(351) = {204, 203};
Line(352) = {203, 201};
Line(353) = {201, 202};
Line(354) = {104, 103};
Line(355) = {103, 101};
Line(356) = {101, 102};
Line(357) = {106, 202};
Line(358) = {105, 201};
Line(359) = {107, 203};
Line(360) = {108, 204};
Line(361) = {204, 104};
Line(362) = {202, 102};
Line(363) = {201, 101};
Line(364) = {203, 103};
Line Loop(306) = {346, -357, -345, 360};
Plane Surface(351) = {306};
Line Loop(307) = {348, 358, 353, -357};
Plane Surface(352) = {307};
Line Loop(308) = {349, 359, 352, -358};
Plane Surface(353) = {308};
Line Loop(309) = {350, 360, 351, -359};
Plane Surface(354) = {309};
Line Loop(310) = {348, 349, 350, 345};
Plane Surface(355) = {310};
Line Loop(311) = {346, -353, -352, -351};
Line Loop(312) = {362, -347, -361, 346};
Plane Surface(357) = {312};
Line Loop(313) = {354, -364, -351, 361};
Plane Surface(358) = {313};
Line Loop(314) = {363, -355, -364, 352};
Plane Surface(359) = {314};
Line Loop(315) = {356, -347, 354, 355};
Plane Surface(361) = {15};
Plane Surface(362) = {17};
Plane Surface(363) = {16};
Plane Surface(364) = {18};
Plane Surface(366) = {14, 310};
Plane Surface(368) = {305, 311};
Line Loop(316) = {332, 333, 330, 331};
Plane Surface(369) = {13, 316};
Line Loop(317) = {311, 308, 309, 310};
Plane Surface(370) = {315, 317};
Line Loop(318) = {363, 356, -362, -353};
Plane Surface(371) = {318};


// electrodes  (in the core)
k=newp;
%ELECTRODES%

Surface Loop(1) = {351, 352, 353, 354, 368, 306, 355};
Volume(3) = {1};
Surface Loop(2) = {357, 371, 359, 370, 358, 368, 327, 315, 319, 323};
Volume(4) = {2};
Surface Loop(3) = {366, 361, 369, 362, 363, 364, 337, 341, 345, 349, 370, 371, 359, 358, 357, 351, 352, 353, 354};
Volume(5) = {3};
Physical Volume("Volcano") = {2, 1};
Physical Volume("Outer") = {5};
Physical Volume("LayerTop") = {3};
Physical Volume("LayerBottom") = {4};
Physical Surface("SurfaceTop") = {366, 355};
Physical Surface("SurfaceBottom") = {369, 350};
Physical Surface("VerticalFaces") = {362, 361, 364, 363};