Mesh.MshFileVersion = 2.2
 
WidthOuterBox=8000.;
ThicknessOuterBox=3500.;
MeshSizeOuterBox=600.;

WidthInnerBox=4500.;
ThicknessInnerBox=1500.;
MeshSizeInnerBox=150.;


// electrodes
numElectrodesX = 11;
numElectrodesY = 11;
SpacingElectrodes = 300.;

MeshSizeAtElectrodes= 20;


// anomaly 1:
Anomaly1Xoffset=600.;
Anomaly1Yoffset=0;

Anomaly1MeshSize=MeshSizeInnerBox/4;
Anomaly1Depth=300.;
Anomaly1Thickness=600.;
Anomaly1LengthX=900.;
Anomaly1LengthY=900.;

// anomaly 2:

Anomaly2Xoffset=-600.;
Anomaly2Yoffset=-600.;

Anomaly2MeshSize=MeshSizeInnerBox/4;
Anomaly2Depth=300.;
Anomaly2Thickness=600.;
Anomaly2LengthX=900.;
Anomaly2LengthY=900.;

// anomaly 3:

Anomaly3Xoffset=-600.;
Anomaly3Yoffset=600.;

Anomaly3MeshSize=MeshSizeInnerBox/4;
Anomaly3Depth=300.;
Anomaly3Thickness=600.;
Anomaly3LengthX=900.;
Anomaly3LengthY=900.;


lenx =(numElectrodesX-1)*SpacingElectrodes;
leny =(numElectrodesY-1)*SpacingElectrodes;


OuterBox0=-WidthOuterBox/2;
InnerBox0=-WidthInnerBox/2;
Anomaly1X0=-Anomaly1LengthX/2+Anomaly1Xoffset;
Anomaly1Y0=-Anomaly1LengthY/2+Anomaly1Yoffset;
Anomaly2X0=-Anomaly2LengthX/2+Anomaly2Xoffset;
Anomaly2Y0=-Anomaly2LengthY/2+Anomaly2Yoffset;
Anomaly3X0=-Anomaly3LengthX/2+Anomaly3Xoffset;
Anomaly3Y0=-Anomaly3LengthY/2+Anomaly3Yoffset;


//outer box points
Point(1) = {-OuterBox0, -OuterBox0, -ThicknessOuterBox, MeshSizeOuterBox};
Point(2) = { OuterBox0, -OuterBox0, -ThicknessOuterBox, MeshSizeOuterBox};
Point(3) = {-OuterBox0,  OuterBox0, -ThicknessOuterBox, MeshSizeOuterBox};
Point(4) = { OuterBox0,  OuterBox0, -ThicknessOuterBox, MeshSizeOuterBox};

Point(5) = {-OuterBox0, -OuterBox0, 0, MeshSizeOuterBox};
Point(6) = { OuterBox0, -OuterBox0, 0, MeshSizeOuterBox};
Point(7) = {-OuterBox0,  OuterBox0, 0, MeshSizeOuterBox};
Point(8) = { OuterBox0,  OuterBox0, 0, MeshSizeOuterBox};

//inner box points
Point(101) = {-InnerBox0, -InnerBox0, -ThicknessInnerBox, MeshSizeInnerBox};
Point(102) = { InnerBox0, -InnerBox0, -ThicknessInnerBox, MeshSizeInnerBox};
Point(103) = {-InnerBox0,  InnerBox0, -ThicknessInnerBox, MeshSizeInnerBox};
Point(104) = { InnerBox0,  InnerBox0, -ThicknessInnerBox, MeshSizeInnerBox};

Point(105) = {-InnerBox0, -InnerBox0, 0, MeshSizeInnerBox};
Point(106) = { InnerBox0, -InnerBox0, 0, MeshSizeInnerBox};
Point(107) = {-InnerBox0,  InnerBox0, 0, MeshSizeInnerBox};
Point(108) = { InnerBox0,  InnerBox0, 0, MeshSizeInnerBox};

//anomaly points
Point(201) = { Anomaly1X0, Anomaly1Y0, -Anomaly1Depth-Anomaly1Thickness, Anomaly1MeshSize};
Point(202) = { Anomaly1X0+Anomaly1LengthX, Anomaly1Y0, -Anomaly1Depth-Anomaly1Thickness, Anomaly1MeshSize};
Point(203) = { Anomaly1X0, Anomaly1Y0+Anomaly1LengthY, -Anomaly1Depth-Anomaly1Thickness, Anomaly1MeshSize};
Point(204) = { Anomaly1X0+Anomaly1LengthX, Anomaly1Y0+Anomaly1LengthY, -Anomaly1Depth-Anomaly1Thickness, Anomaly1MeshSize};

Point(205) = { Anomaly1X0, Anomaly1Y0, -Anomaly1Depth, Anomaly1MeshSize};
Point(206) = { Anomaly1X0+Anomaly1LengthX, Anomaly1Y0, -Anomaly1Depth, Anomaly1MeshSize};
Point(207) = { Anomaly1X0, Anomaly1Y0+Anomaly1LengthY, -Anomaly1Depth, Anomaly1MeshSize};
Point(208) = { Anomaly1X0+Anomaly1LengthX, Anomaly1Y0+Anomaly1LengthY, -Anomaly1Depth, Anomaly1MeshSize};


Point(301) = { Anomaly2X0, Anomaly2Y0, -Anomaly2Depth-Anomaly2Thickness, Anomaly2MeshSize};
Point(302) = { Anomaly2X0+Anomaly2LengthX, Anomaly2Y0, -Anomaly2Depth-Anomaly2Thickness, Anomaly2MeshSize};
Point(303) = { Anomaly2X0, Anomaly2Y0+Anomaly2LengthY, -Anomaly2Depth-Anomaly2Thickness, Anomaly2MeshSize};
Point(304) = { Anomaly2X0+Anomaly2LengthX, Anomaly2Y0+Anomaly2LengthY, -Anomaly2Depth-Anomaly2Thickness, Anomaly2MeshSize};

Point(305) = { Anomaly2X0, Anomaly2Y0, -Anomaly2Depth, Anomaly2MeshSize};
Point(306) = { Anomaly2X0+Anomaly2LengthX, Anomaly2Y0, -Anomaly2Depth, Anomaly2MeshSize};
Point(307) = { Anomaly2X0, Anomaly2Y0+Anomaly2LengthY, -Anomaly2Depth, Anomaly2MeshSize};
Point(308) = { Anomaly2X0+Anomaly2LengthX, Anomaly2Y0+Anomaly2LengthY, -Anomaly2Depth, Anomaly2MeshSize};

Point(401) = { Anomaly3X0, Anomaly3Y0, -Anomaly3Depth-Anomaly3Thickness, Anomaly3MeshSize};
Point(402) = { Anomaly3X0+Anomaly3LengthX, Anomaly3Y0, -Anomaly3Depth-Anomaly3Thickness, Anomaly3MeshSize};
Point(403) = { Anomaly3X0, Anomaly3Y0+Anomaly3LengthY, -Anomaly3Depth-Anomaly3Thickness, Anomaly3MeshSize};
Point(404) = { Anomaly3X0+Anomaly3LengthX, Anomaly3Y0+Anomaly3LengthY, -Anomaly3Depth-Anomaly3Thickness, Anomaly3MeshSize};

Point(405) = { Anomaly3X0, Anomaly3Y0, -Anomaly3Depth, Anomaly3MeshSize};
Point(406) = { Anomaly3X0+Anomaly3LengthX, Anomaly3Y0, -Anomaly3Depth, Anomaly3MeshSize};
Point(407) = { Anomaly3X0, Anomaly3Y0+Anomaly3LengthY, -Anomaly3Depth, Anomaly3MeshSize};
Point(408) = { Anomaly3X0+Anomaly3LengthX, Anomaly3Y0+Anomaly3LengthY, -Anomaly3Depth, Anomaly3MeshSize};

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


//make inner box lines
//bottom
Line(101) = {101,102} ;
Line(102) = {103,104} ;
Line(103) = {101,103} ;
Line(104) = {102,104} ;
//top
Line(105) = {105,106} ;
Line(106) = {107,108} ;
Line(107) = {105,107} ;
Line(108) = {106,108} ;
//sides
Line(109) = {101,105} ;
Line(110) = {102,106} ;
Line(111) = {103,107} ;
Line(112) = {104,108} ;

// outward normals!!! bottom top front back left right
Line Loop(113) = {103,102,-104,-101} ;   
Line Loop(114) = {105,108,-106,-107} ;   
Line Loop(115) = {101,110,-105,-109} ;  
Line Loop(116) = {111,106,-112,-102} ; 
Line Loop(117) = {109,107,-111,-103} ;  
Line Loop(118) = {104,112,-108,-110} ; 

//make anomaly lines1 
//bottom
Line(201) = {201,202} ;
Line(202) = {203,204} ;
Line(203) = {201,203} ;
Line(204) = {202,204} ;
//top
Line(205) = {205,206} ;
Line(206) = {207,208} ;
Line(207) = {205,207} ;
Line(208) = {206,208} ;
//sides
Line(209) = {201,205} ;
Line(210) = {202,206} ;
Line(211) = {203,207} ;
Line(212) = {204,208} ;

// outward normals!!! bottom top front back left right
Line Loop(213) = {203,202,-204,-201} ;   
Line Loop(214) = {205,208,-206,-207} ;   
Line Loop(215) = {201,210,-205,-209} ;  
Line Loop(216) = {211,206,-212,-202} ; 
Line Loop(217) = {209,207,-211,-203} ;  
Line Loop(218) = {204,212,-208,-210} ; 

//make anomaly lines 2
//bottom
Line(301) = {301,302} ;
Line(302) = {303,304} ;
Line(303) = {301,303} ;
Line(304) = {302,304} ;
//top
Line(305) = {305,306} ;
Line(306) = {307,308} ;
Line(307) = {305,307} ;
Line(308) = {306,308} ;
//sides
Line(309) = {301,305} ;
Line(310) = {302,306} ;
Line(311) = {303,307} ;
Line(312) = {304,308} ;

// outward normals!!! bottom top front back left right
Line Loop(313) = {303,302,-304,-301} ;   
Line Loop(314) = {305,308,-306,-307} ;   
Line Loop(315) = {301,310,-305,-309} ;  
Line Loop(316) = {311,306,-312,-302} ; 
Line Loop(317) = {309,307,-311,-303} ;  
Line Loop(318) = {304,312,-308,-310} ; 


//make anomaly lines 3
//bottom
Line(401) = {401,402} ;
Line(402) = {403,404} ;
Line(403) = {401,403} ;
Line(404) = {402,404} ;
//top
Line(405) = {405,406} ;
Line(406) = {407,408} ;
Line(407) = {405,407} ;
Line(408) = {406,408} ;
//sides
Line(409) = {401,405} ;
Line(410) = {402,406} ;
Line(411) = {403,407} ;
Line(412) = {404,408} ;

// outward normals!!! bottom top front back left right
Line Loop(413) = {403,402,-404,-401} ;   
Line Loop(414) = {405,408,-406,-407} ;   
Line Loop(415) = {401,410,-405,-409} ;  
Line Loop(416) = {411,406,-412,-402} ; 
Line Loop(417) = {409,407,-411,-403} ;  
Line Loop(418) = {404,412,-408,-410} ; 

// Surfaces
Plane Surface(101) = {113};  
Plane Surface(102) = {114};  
Plane Surface(103) = {115};  
Plane Surface(104) = {116};  
Plane Surface(105) = {117};  
Plane Surface(106) = {118};  

Plane Surface(1) = {13};  
Plane Surface(2) = {14, 114};  
Plane Surface(3) = {15};  
Plane Surface(4) = {16};  
Plane Surface(5) = {17};  
Plane Surface(6) = {18};  


Plane Surface(201) = {213};  
Plane Surface(202) = {214};  
Plane Surface(203) = {215};  
Plane Surface(204) = {216};  
Plane Surface(205) = {217};  
Plane Surface(206) = {218}; 

Plane Surface(301) = {313};  
Plane Surface(302) = {314};  
Plane Surface(303) = {315};  
Plane Surface(304) = {316};  
Plane Surface(305) = {317};  
Plane Surface(306) = {318}; 

Plane Surface(401) = {413};  
Plane Surface(402) = {414};  
Plane Surface(403) = {415};  
Plane Surface(404) = {416};  
Plane Surface(405) = {417};  
Plane Surface(406) = {418}; 

//electrodes
k=newp;
Point(k+1)={ ( (-1.)/(numElectrodesX-1)-0.5)*lenx, 0.0, 0.0, MeshSizeAtElectrodes};
Point{k+1} In Surface{102};
Point(k+2)={  ( (numElectrodesX+0.)/(numElectrodesX-1)-0.5)*lenx, 0., 0.0, MeshSizeAtElectrodes};
Point{k+2} In Surface{102};
Point(k+3)={ 0., ((-1.)/(numElectrodesY-1) -0.5)*leny, 0.0, MeshSizeAtElectrodes};
Point{k+3} In Surface{102};
Point(k+4)={ 0.,  ((numElectrodesY+0.)/(numElectrodesY-1) -0.5)*leny, 0.0, MeshSizeAtElectrodes};
Point{k+4} In Surface{102};


k=newp;
For i In{0:numElectrodesY-1}
  For j In{0:numElectrodesX-1}
    xpos=(j/(numElectrodesX-1) -0.5)*lenx;
    ypos=(i/(numElectrodesY-1) -0.5)*leny;
    Point(k+i*numElectrodesX+j)={ xpos, ypos, 0.0, MeshSizeAtElectrodes};
    Point{k+i*numElectrodesX+j} In Surface{102};
  EndFor
EndFor
Coherence;


Surface Loop(220) = {201,202,203,204,205,206};
Volume(20) = {220};
Physical Volume("Anomaly1")={20};

Surface Loop(320) = {301,302,303,304,305,306};
Volume(30) = {320};
Physical Volume("Anomaly2")={30};

Surface Loop(420) = {401,402,403,404,405,406};
Volume(40) = {420};
Physical Volume("Anomaly3")={40};


Surface Loop(120) = {101,102,103,104,105,106};
Volume(10) = {120, 220, 320, 420};

Physical Volume("InnerBox")={10};

Surface Loop(20) = {1,2,3,4,5,6, -101,-103,-104,-105,-106};
Volume(1) = {20};
Physical Volume("OuterBox")={1};

