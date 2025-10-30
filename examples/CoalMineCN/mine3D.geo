Mesh.MshFileVersion = 2.2;
// ... geometry:
GalleryHeight = 5;
SeamHeight = 13;
CoreHeight = 100;
BaseDepth = 100;
PaddingX=1500;
PaddingY=1500;
PaddingZ=1500;
//... center coordinates in original IGES data
OriginalCenterX = 496192.29259967303;
OriginalCenterY = 890424.0922468755;
//... mesh sizes ...
meshSizeRemoteCore = 15.;
meshSizeCore = 5.;
meshSizeCloseToSeam = 3;
meshSizeSeam = 2.5/1.2;
meshSizeElectrodes = 2.5/400;
meshSizePadding = 2*(241.06384326395346+PaddingX)/15;
//... mine lay out
Point(226) = {-44.35527517402079, 56.828797810478136, 0., meshSizeSeam};
Point(312) = {147.4847248259466 + 4., 56.82851888751611, 0., meshSizeSeam};
Point(314) = {147.48472712998046 +4., -143.17124071449507, 0., meshSizeSeam};
Point(318) = {-65.57938983303029, -178.42203949647956, 0., meshSizeCloseToSeam};
Point(319) = {-65.57938983303029, -163.74508660647552, 0., meshSizeSeam};
Point(320) = {-60.57938983303029, -178.42203949647956, 0., meshSizeSeam};
Point(321) = {-60.57938983303029, -163.74508660647552, 0., meshSizeSeam};
Point(400) = {69.35798970394535, -148.57202080648858, 0., meshSizeCloseToSeam};
Point(401) = {64.35798970394535, -148.57174933550414, 0., meshSizeCloseToSeam};
Point(403) = {197.48472482495708, 56.82837679947261, 0., meshSizeCore};
Point(404) = {197.48471409198828, -143.17125238850713, 0., meshSizeCore};
Point(407) = {-20.043531734030694, -143.17124071449507, 0., meshSizeSeam};
Point(411) = {-60.4151782640256, 72.22858866746537, 0., meshSizeSeam};
Point(415) = {-115.23431354801869, 77.62831308343448, 0., meshSizeSeam};
Point(420) = {-115.23432427900843, -183.82199773250613, 0., meshSizeCore};
Point(421) = {-115.23432427900843, -178.4219617624767, 0., meshSizeCloseToSeam};
Point(426) = {-115.23432427900843, -148.57135686755646, 0., meshSizeCloseToSeam};
Point(427) = {-115.23432427900843, -143.17110019456595, 0., meshSizeSeam};
Point(428) = {-120.63432427705266, -148.57135686755646, 0., meshSizeCloseToSeam};
Point(432) = {-115.23432427900843, 72.22831308352761, 0., meshSizeSeam};
Point(434) = {204.35999706096482, -178.5718327035429, 0., meshSizeCore};
Point(436) = {69.35798970394535, -178.5718327035429, 0., meshSizeCore};
Point(438) = {64.35798970394535, -178.5718327035429, 0., meshSizeCore};
Point(440) = {209.35999706096482, -183.8216939845588, 0., meshSizeCore};
Point(444) = {209.35999706096482, -148.5718042094959, 0., meshSizeCore};
Point(445) = {241.06384326395346, -148.57180284056813, 0., meshSizeCore};
Point(446) = {-47.65758156304946, -178.5718327035429, 0., meshSizeCloseToSeam};
Point(450) = {204.35999706096482, -148.57214718847536, 0., meshSizeCloseToSeam};
Point(451) = {-17.80705362104345, -148.57151584257372, 0., meshSizeCloseToSeam};
Point(455) = {-55.29433802800486, -178.4220470135333, 0., meshSizeSeam};
Point(459) = {-194.93432427401422, -148.5713807165157, 0., meshSizeCore};
Point(460) = {-25.443807955016382, -148.57151693548076, 0., meshSizeSeam};
Point(461) = {-194.93432427401422, -143.17123132850975, 0., meshSizeCore};
Point(462) = {241.06384326395346, -143.17125472647604, 0., meshSizeCore};
Point(463) = {-42.77832962805405, 62.22849326150026, 0., meshSizeCloseToSeam};
Point(464) = {241.06384326395346, 62.22809668444097, 0., meshSizeCore};
Point(465) = {241.06384326395346, 56.82842257351149, 0., meshSizeCore};
Point(466) = {-194.93432427401422, 56.828823404503055, 0., meshSizeCore};
Point(467) = {-230.0899094590568, -233.97, 0., meshSizeCore};
Point(468) = {-230.0899094590568, 112.3783179235179, 0., meshSizeCore};
Point(471) = {-120.63432427705266, 62.22865014243871, 0., meshSizeCloseToSeam};
Point(472) = {-120.63432427705266, 97.22832381352782, 0., meshSizeCore};
Point(473) = {-120.63432427705266, -33.97255319147371, 0., meshSizeCloseToSeam};
Point(474) = {-120.63432427705266, 56.82856875646394, 0., meshSizeCloseToSeam};
Point(475) = {-120.63432427705266, 36.828050636453554, 0., meshSizeCloseToSeam};
Point(476) = {-120.63432427705266, 42.228050637524575, 0., meshSizeCloseToSeam};
Point(477) = {-115.23432427900843, 102.62831308343448, 0., meshSizeCore};
Point(478) = {-115.23432427900843, 112.3783179235179, 0., meshSizeCore};
Point(479) = {-115.23432427900843, 97.22847146447748, 0., meshSizeCore};
Point(484) = {-120.63432427705266, 112.3783179235179, 0., meshSizeCore};
Point(485) = {-120.63432427705266, 102.62832381448243, 0., meshSizeCore};
Point(486) = {-120.63432427705266, -143.17110907949973, 0., meshSizeCloseToSeam};
Point(487) = {-120.63432427705266, -39.37255319149699, 0., meshSizeCloseToSeam};
Point(488) = {-120.63432427705266, -233.97, 0., meshSizeCore};
Point(489) = {-120.63432427705266, -168.57194936356973, 0., meshSizeCore};
Point(490) = {-120.63432427705266, -163.17194929055404, 0., meshSizeCore};
Point(497) = {-51.691298495046794, 77.62818681052886, 0., meshSizeCore};
Point(498) = {-32.09133736306103, 97.22814785945229, 0., meshSizeCore};
Point(500) = {241.06384326395346, 102.62759396550246, 0., meshSizeCore};
Point(501) = {241.06384326395346, 97.22847994044423, 0., meshSizeCore};
Point(502) = {130.76615817594575, 97.22860790649429, 0., meshSizeCore};
Point(504) = {241.06384326395346, 72.6776050424669, 0., meshSizeCore};
Point(506) = {110.71541434299434, 77.17786407051608, 0., meshSizeCore};
Point(507) = {241.06384326395346, 77.17768545344006, 0., meshSizeCore};
Point(512) = {-194.93408996105427, 102.62847146543209, 0., meshSizeCore};
Point(513) = {32.585906428983435, 77.62773523549549, 0., meshSizeCore};
Point(514) = {103.52839129196946, 77.62759425945114, 0., meshSizeCore};
Point(515) = {32.58578806096921, 97.22801933146548, 0., meshSizeCore};
Point(517) = {27.185534243937582, 97.22803006449249, 0., meshSizeCore};
Point(518) = {27.185906429949682, 77.62847146450076, 0., meshSizeCore};
Point(529) = {-52.778399760019965, 72.22856339346617, 0., meshSizeCore};
Point(530) = {-44.054221655009314, 77.62851041345857, 0., meshSizeCore};
Point(531) = {-24.45426052302355, 97.22847146447748, 0., meshSizeCore};
Point(532) = {-200.9352204600582, 112.3783179235179, 0., meshSizeCore};
Point(533) = {-200.93421152001247, -233.97, 0., meshSizeCore};
Point(535) = {-160.5343242760282, -143.17117472854443, 0., meshSizeCore};
Point(536) = {-160.5343242760282, -38.37255319149699, 0., meshSizeCore};
Point(537) = {-160.5343242760282, -233.97, 0., meshSizeCore};
Point(538) = {-160.5343242760282, -148.57128977356479, 0., meshSizeCore};
Point(539) = {-155.33432427304797, -143.17116617353167, 0., meshSizeCore};
Point(540) = {-155.33432427304797, -39.37255319149699, 0., meshSizeCore};
Point(541) = {-155.33432427304797, -233.97, 0., meshSizeCore};
Point(542) = {-155.33432427304797, -148.57129851856735, 0., meshSizeCore};
Point(543) = {-194.93432427401422, -38.37255319149699, 0., meshSizeCore};
Point(544) = {-160.5343242760282, -34.97255319147371, 0., meshSizeCore};
Point(545) = {-194.93432427401422, -34.97255319147371, 0., meshSizeCore};
Point(546) = {-155.33432427304797, -33.97255319147371, 0., meshSizeCore};
Point(548) = {-194.93432427401422, 97.22847146447748, 0., meshSizeCore};
Point(550) = {123.12915126694134, 97.22835423343349, 0., meshSizeCore};
Point(559) = {-194.93432427401422, -233.97, 0., meshSizeCore};
Point(561) = {-194.93432427401422, 62.2287680885056, 0., meshSizeCore};
Point(562) = {-194.93432427401422, 112.3783179235179, 0., meshSizeCore};
Point(564) = {-160.5343242760282, 56.82867569348309, 0., meshSizeCore};
Point(565) = {-160.5343242760282, 97.22840310353786, 0., meshSizeCore};
Point(566) = {-160.5343242760282, 62.22871348145418, 0., meshSizeCore};
Point(567) = {-160.5343242760282, 112.3783179235179, 0., meshSizeCore};
Point(568) = {-160.5343242760282, 102.62840310449246, 0., meshSizeCore};
Point(569) = {-155.33432427304797, 56.82866175752133, 0., meshSizeCore};
Point(570) = {-155.33432427304797, 97.22839277051389, 0., meshSizeCore};
Point(571) = {-155.33432427304797, 62.22870522644371, 0., meshSizeCore};
Point(572) = {-155.33432427304797, 112.3783179235179, 0., meshSizeCore};
Point(573) = {-155.33432427304797, 102.62839277053718, 0., meshSizeCore};
Point(574) = {-194.93432427401422, -168.5718436865136, 0., meshSizeCore};
Point(575) = {-194.93432427401422, -163.171843685559, 0., meshSizeCore};
Point(576) = {-194.93432427401422, 36.82815631350968, 0., meshSizeCore};
Point(577) = {-194.93432427401422, 42.228156314464286, 0., meshSizeCore};
Point(578) = {-155.33432427304797, -168.57190000952687, 0., meshSizeCore};
Point(579) = {-155.33432427304797, -163.17189997050446, 0., meshSizeCore};
Point(580) = {-155.33432427304797, 36.828099991427734, 0., meshSizeCore};
Point(581) = {-155.33432427304797, 42.22809999145102, 0., meshSizeCore};
Point(583) = {-160.5343242760282, -168.57189261447638, 0., meshSizeCore};
Point(584) = {-160.5343242760282, -163.1718925795285, 0., meshSizeCore};
Point(585) = {-160.5343242760282, 36.82810738647822, 0., meshSizeCore};
Point(586) = {-160.5343242760282, 42.22810738743283, 0., meshSizeCore};
Point(588) = {241.06384326395346, 112.378317923517, 0, meshSizeCore};
Point(589) = {241.06384326395346, -233.97, 0, meshSizeCore};
Point(590) = {-115.23432427900843, -233.97, 0, meshSizeCore};
Point(591) = {147.4847248259466, 62.22849326150026, 0., meshSizeCloseToSeam};
Point(592) = {197.4847140919882, 62.22849326150026, 0., meshSizeCore};
Point(593) = {147.48472712998046, -148.5720208064885, 0., meshSizeCloseToSeam};
Point(594) = {197.4847140919882, -148.5720208064885, 0., meshSizeCore};

Line(1) = {532, 533};
Line(2) = {562, 512};
Line(3) = {512, 568};
Line(4) = {568, 567};
Line(5) = {548, 565};
Line(6) = {565, 566};
Line(7) = {566, 561};
Line(8) = {561, 548};
Line(9) = {466, 564};
Line(10) = {564, 586};
Line(11) = {586, 577};
Line(12) = {577, 466};
Line(13) = {585, 576};
Line(14) = {585, 544};
Line(15) = {544, 545};
Line(16) = {545, 576};
Line(17) = {543, 536};
Line(18) = {536, 535};
Line(19) = {535, 461};
Line(20) = {461, 543};
Line(21) = {459, 538};
Line(22) = {538, 584};
Line(23) = {584, 575};
Line(25) = {575, 459};
Line(26) = {574, 583};
Line(27) = {583, 537};
Line(28) = {574, 559};
Line(29) = {572, 573};
Line(30) = {573, 485};
Line(31) = {485, 484};
Line(32) = {570, 472};
Line(33) = {472, 471};
Line(34) = {471, 571};
Line(35) = {571, 570};
Line(36) = {569, 474};
Line(37) = {474, 476};
Line(38) = {476, 581};
Line(39) = {569, 581};
Line(40) = {580, 475};
Line(41) = {475, 473};
Line(42) = {473, 546};
Line(43) = {546, 580};
Line(44) = {540, 487};
Line(45) = {487, 486};
Line(46) = {486, 539};
Line(47) = {539, 540};
Line(48) = {542, 428};
Line(49) = {428, 490};
Line(50) = {490, 579};
Line(51) = {579, 542};
Line(52) = {578, 489};
Line(53) = {489, 488};
Line(54) = {488, 541};
Line(55) = {578, 541};
Line(56) = {537, 541};
Line(57) = {537, 559};
Line(58) = {559, 533};
Line(59) = {533, 467};
Line(60) = {467, 468};
Line(61) = {468, 532};
Line(62) = {532, 562};
Line(63) = {562, 567};
Line(64) = {567, 572};
Line(65) = {572, 484};
Line(66) = {484, 478};
Line(67) = {478, 477};
Line(69) = {477, 500};
Line(70) = {479, 498};
Line(71) = {531, 517};
Line(72) = {515, 550};
Line(73) = {502, 501};
Line(74) = {501, 500};
Line(75) = {502, 506};
Line(76) = {506, 507};
Line(77) = {550, 514};
Line(78) = {514, 513};
Line(79) = {513, 515};
Line(80) = {517, 518};
Line(81) = {518, 530};
Line(82) = {530, 531};
Line(83) = {497, 498};
Line(84) = {415, 497};
Line(87) = {507, 504};
Line(88) = {504, 529};
Line(89) = {529, 463};
Line(91) = {464, 504};
Line(92) = {507, 501};
Line(93) = {226, 411};
Line(94) = {411, 432};
Line(95) = {415, 479};
Line(96) = {432, 427};
Line(99) = {465, 464};
Line(100) = {312, 314};
Line(101) = {465, 462};
Line(102) = {403, 404};
Line(103) = {404, 462};
Line(104) = {445, 444};
Line(105) = {444, 440};
Line(106) = {440, 420};
Line(109) = {407, 314};
Line(110) = {455, 460};
Line(111) = {451, 446};
Line(112) = {427, 407};
Line(113) = {426, 460};
Line(115) = {434, 450};
Line(116) = {426, 421};
Line(117) = {401, 438};
Line(118) = {400, 436};
Line(120) = {436, 434};
Line(121) = {438, 446};
Line(122) = {451, 401};
Line(123) = {421, 318};
Line(124) = {318, 319};
Line(125) = {319, 321};
Line(126) = {321, 320};
Line(127) = {320, 455};
Line(128) = {226, 312};
Line(129) = {403, 465};
Line(130) = {462, 445};
Line(131) = {488, 590};
Line(132) = {590, 589};
Line(133) = {589, 445};
Line(134) = {420, 590};
Line(135) = {478, 588};
Line(136) = {500, 588};
Line(137) = {463, 591};
Line(138) = {591, 592};
Line(139) = {592, 464};
Line(140) = {400, 593};
Line(141) = {593, 594};
Line(142) = {594, 450};



Curve Loop(1) = {60, 61, 1, 59};
Plane Surface(1) = {1};
Curve Loop(2) = {63, -4, -3, -2};
Plane Surface(2) = {2};
Curve Loop(3) = {5, 6, 7, 8};
Plane Surface(3) = {3};
Curve Loop(4) = {9, 10, 11, 12};
Plane Surface(4) = {4};
Curve Loop(5) = {13, -16, -15, -14};
Plane Surface(5) = {5};
Curve Loop(6) = {17, 18, 19, 20};
Plane Surface(6) = {6};
Curve Loop(7) = {21, 22, 23, 25};
Plane Surface(7) = {7};
Curve Loop(8) = {26, 27, 57, -28};
Plane Surface(8) = {8};
Curve Loop(9) = {52, 53, 54, -55};
Plane Surface(9) = {9};
Curve Loop(10) = {48, 49, 50, 51};
Plane Surface(10) = {10};
Curve Loop(11) = {44, 45, 46, 47};
Plane Surface(11) = {11};
Curve Loop(12) = {40, 41, 42, 43};
Plane Surface(12) = {12};
Curve Loop(13) = {36, 37, 38, -39};
Plane Surface(13) = {13};
Curve Loop(14) = {32, 33, 34, 35};
Plane Surface(14) = {14};
Curve Loop(15) = {65, -31, -30, -29};
Plane Surface(15) = {15};
Curve Loop(16) = {135, -136, -69, -67};
Plane Surface(16) = {16};
Curve Loop(17) = {70, -83, -84, 95};
Plane Surface(17) = {17};
Curve Loop(18) = {71, 80, 81, 82};
Plane Surface(18) = {18};
Curve Loop(19) = {72, 77, 78, 79};
Plane Surface(19) = {19};
Curve Loop(20) = {73, -92, -76, -75};
Plane Surface(20) = {20};
Curve Loop(21) = {94, 96, 112, 109, -100, -128, 93};
Plane Surface(21) = {21};
Curve Loop(23) = {101, -103, -102, 129};
Plane Surface(23) = {23};
Curve Loop(24) = {104, 105, 106, 134, 132, 133};
Plane Surface(24) = {24};
Curve Loop(26) = {122, 117, 121, -111};
Plane Surface(26) = {26};
Curve Loop(27) = {113, -110, -127, -126, -125, -124, -123, -116};
Plane Surface(27) = {27};
Curve Loop(28) = {137, 138, 139, 91, 88, 89};
Plane Surface(28) = {28};
Curve Loop(29) = {120, 115, -142, -141, -140, 118};
Plane Surface(29) = {29};
Curve Loop(30) = {1, -58, -28, 26, 27, 56, -55, 52, 53, 131, -134, -106, -105, -104, -130, -103, -102, 129, 99, -139, -138, -137, -89, -88, -87, -76, -75, 73, 74, -69, -67, -66, -31, -30, -29, -64, -4, -3, -2, -62};
Plane Surface(30) = {3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 17, 18, 19, 21, 26, 27, 29, 30};
// ... top of galleries ... 
Point(821) = {-44.35527517402079, 56.828797810478136, GalleryHeight, meshSizeSeam};
Point(907) = {147.4847248259466 + 4., 56.82851888751611, GalleryHeight, meshSizeSeam};
Point(909) = {147.48472712998046 + 4., -143.17124071449507, GalleryHeight, meshSizeSeam};
Point(913) = {-65.57938983303029, -178.42203949647956, GalleryHeight, meshSizeCloseToSeam};
Point(914) = {-65.57938983303029, -163.74508660647552, GalleryHeight, meshSizeSeam};
Point(915) = {-60.57938983303029, -178.42203949647956, GalleryHeight, meshSizeSeam};
Point(916) = {-60.57938983303029, -163.74508660647552, GalleryHeight, meshSizeSeam};
Point(995) = {69.35798970394535, -148.57202080648858, GalleryHeight, meshSizeCloseToSeam};
Point(996) = {64.35798970394535, -148.57174933550414, GalleryHeight, meshSizeCloseToSeam};
Point(998) = {197.48472482495708, 56.82837679947261, GalleryHeight, meshSizeCore};
Point(999) = {197.48471409198828, -143.17125238850713, GalleryHeight, meshSizeCore};
Point(1002) = {-20.043531734030694, -143.17124071449507, GalleryHeight, meshSizeSeam};
Point(1006) = {-60.4151782640256, 72.22858866746537, GalleryHeight, meshSizeSeam};
Point(1010) = {-115.23431354801869, 77.62831308343448, GalleryHeight, meshSizeSeam};
Point(1015) = {-115.23432427900843, -183.82199773250613, GalleryHeight, meshSizeCore};
Point(1016) = {-115.23432427900843, -178.4219617624767, GalleryHeight, meshSizeCloseToSeam};
Point(1021) = {-115.23432427900843, -148.57135686755646, GalleryHeight, meshSizeCloseToSeam};
Point(1022) = {-115.23432427900843, -143.17110019456595, GalleryHeight, meshSizeSeam};
Point(1023) = {-120.63432427705266, -148.57135686755646, GalleryHeight, meshSizeCloseToSeam};
Point(1027) = {-115.23432427900843, 72.22831308352761, GalleryHeight, meshSizeSeam};
Point(1029) = {204.35999706096482, -178.5718327035429, GalleryHeight, meshSizeCore};
Point(1031) = {69.35798970394535, -178.5718327035429, GalleryHeight, meshSizeCore};
Point(1033) = {64.35798970394535, -178.5718327035429, GalleryHeight, meshSizeCore};
Point(1035) = {209.35999706096482, -183.8216939845588, GalleryHeight, meshSizeCore};
Point(1039) = {209.35999706096482, -148.5718042094959, GalleryHeight, meshSizeCore};
Point(1040) = {241.06384326395346, -148.57180284056813, GalleryHeight, meshSizeCore};
Point(1041) = {-47.65758156304946, -178.5718327035429, GalleryHeight, meshSizeCloseToSeam};
Point(1045) = {204.35999706096482, -148.57214718847536, GalleryHeight, meshSizeCloseToSeam};
Point(1046) = {-17.80705362104345, -148.57151584257372, GalleryHeight, meshSizeCloseToSeam};
Point(1050) = {-55.29433802800486, -178.4220470135333, GalleryHeight, meshSizeSeam};
Point(1054) = {-194.93432427401422, -148.5713807165157, GalleryHeight, meshSizeCore};
Point(1055) = {-25.443807955016382, -148.57151693548076, GalleryHeight, meshSizeSeam};
Point(1056) = {-194.93432427401422, -143.17123132850975, GalleryHeight, meshSizeCore};
Point(1057) = {241.06384326395346, -143.17125472647604, GalleryHeight, meshSizeCore};
Point(1058) = {-42.77832962805405, 62.22849326150026, GalleryHeight, meshSizeCloseToSeam};
Point(1059) = {241.06384326395346, 62.22809668444097, GalleryHeight, meshSizeCore};
Point(1060) = {241.06384326395346, 56.82842257351149, GalleryHeight, meshSizeCore};
Point(1061) = {-194.93432427401422, 56.828823404503055, GalleryHeight, meshSizeCore};
Point(1066) = {-120.63432427705266, 62.22865014243871, GalleryHeight, meshSizeCloseToSeam};
Point(1067) = {-120.63432427705266, 97.22832381352782, GalleryHeight, meshSizeCore};
Point(1068) = {-120.63432427705266, -33.97255319147371, GalleryHeight, meshSizeCloseToSeam};
Point(1069) = {-120.63432427705266, 56.82856875646394, GalleryHeight, meshSizeCloseToSeam};
Point(1070) = {-120.63432427705266, 36.828050636453554, GalleryHeight, meshSizeCloseToSeam};
Point(1071) = {-120.63432427705266, 42.228050637524575, GalleryHeight, meshSizeCloseToSeam};
Point(1072) = {-115.23432427900843, 102.62831308343448, GalleryHeight, meshSizeCore};
Point(1073) = {-115.23432427900843, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1074) = {-115.23432427900843, 97.22847146447748, GalleryHeight, meshSizeCore};
Point(1079) = {-120.63432427705266, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1080) = {-120.63432427705266, 102.62832381448243, GalleryHeight, meshSizeCore};
Point(1081) = {-120.63432427705266, -143.17110907949973, GalleryHeight, meshSizeCloseToSeam};
Point(1082) = {-120.63432427705266, -39.37255319149699, GalleryHeight, meshSizeCloseToSeam};
Point(1083) = {-120.63432427705266, -233.97, GalleryHeight, meshSizeCore};
Point(1084) = {-120.63432427705266, -168.57194936356973, GalleryHeight, meshSizeCore};
Point(1085) = {-120.63432427705266, -163.17194929055404, GalleryHeight, meshSizeCore};
Point(1092) = {-51.691298495046794, 77.62818681052886, GalleryHeight, meshSizeCore};
Point(1093) = {-32.09133736306103, 97.22814785945229, GalleryHeight, meshSizeCore};
Point(1095) = {241.06384326395346, 102.62759396550246, GalleryHeight, meshSizeCore};
Point(1096) = {241.06384326395346, 97.22847994044423, GalleryHeight, meshSizeCore};
Point(1097) = {130.76615817594575, 97.22860790649429, GalleryHeight, meshSizeCore};
Point(1099) = {241.06384326395346, 72.6776050424669, GalleryHeight, meshSizeCore};
Point(1101) = {110.71541434299434, 77.17786407051608, GalleryHeight, meshSizeCore};
Point(1102) = {241.06384326395346, 77.17768545344006, GalleryHeight, meshSizeCore};
Point(1107) = {-194.93408996105427, 102.62847146543209, GalleryHeight, meshSizeCore};
Point(1108) = {32.585906428983435, 77.62773523549549, GalleryHeight, meshSizeCore};
Point(1109) = {103.52839129196946, 77.62759425945114, GalleryHeight, meshSizeCore};
Point(1110) = {32.58578806096921, 97.22801933146548, GalleryHeight, meshSizeCore};
Point(1112) = {27.185534243937582, 97.22803006449249, GalleryHeight, meshSizeCore};
Point(1113) = {27.185906429949682, 77.62847146450076, GalleryHeight, meshSizeCore};
Point(1124) = {-52.778399760019965, 72.22856339346617, GalleryHeight, meshSizeCore};
Point(1125) = {-44.054221655009314, 77.62851041345857, GalleryHeight, meshSizeCore};
Point(1126) = {-24.45426052302355, 97.22847146447748, GalleryHeight, meshSizeCore};
Point(1127) = {-200.9352204600582, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1128) = {-200.93421152001247, -233.97, GalleryHeight, meshSizeCore};
Point(1130) = {-160.5343242760282, -143.17117472854443, GalleryHeight, meshSizeCore};
Point(1131) = {-160.5343242760282, -38.37255319149699, GalleryHeight, meshSizeCore};
Point(1132) = {-160.5343242760282, -233.97, GalleryHeight, meshSizeCore};
Point(1133) = {-160.5343242760282, -148.57128977356479, GalleryHeight, meshSizeCore};
Point(1134) = {-155.33432427304797, -143.17116617353167, GalleryHeight, meshSizeCore};
Point(1135) = {-155.33432427304797, -39.37255319149699, GalleryHeight, meshSizeCore};
Point(1136) = {-155.33432427304797, -233.97, GalleryHeight, meshSizeCore};
Point(1137) = {-155.33432427304797, -148.57129851856735, GalleryHeight, meshSizeCore};
Point(1138) = {-194.93432427401422, -38.37255319149699, GalleryHeight, meshSizeCore};
Point(1139) = {-160.5343242760282, -34.97255319147371, GalleryHeight, meshSizeCore};
Point(1140) = {-194.93432427401422, -34.97255319147371, GalleryHeight, meshSizeCore};
Point(1141) = {-155.33432427304797, -33.97255319147371, GalleryHeight, meshSizeCore};
Point(1143) = {-194.93432427401422, 97.22847146447748, GalleryHeight, meshSizeCore};
Point(1145) = {123.12915126694134, 97.22835423343349, GalleryHeight, meshSizeCore};
Point(1154) = {-194.93432427401422, -233.97, GalleryHeight, meshSizeCore};
Point(1156) = {-194.93432427401422, 62.2287680885056, GalleryHeight, meshSizeCore};
Point(1157) = {-194.93432427401422, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1159) = {-160.5343242760282, 56.82867569348309, GalleryHeight, meshSizeCore};
Point(1160) = {-160.5343242760282, 97.22840310353786, GalleryHeight, meshSizeCore};
Point(1161) = {-160.5343242760282, 62.22871348145418, GalleryHeight, meshSizeCore};
Point(1162) = {-160.5343242760282, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1163) = {-160.5343242760282, 102.62840310449246, GalleryHeight, meshSizeCore};
Point(1164) = {-155.33432427304797, 56.82866175752133, GalleryHeight, meshSizeCore};
Point(1165) = {-155.33432427304797, 97.22839277051389, GalleryHeight, meshSizeCore};
Point(1166) = {-155.33432427304797, 62.22870522644371, GalleryHeight, meshSizeCore};
Point(1167) = {-155.33432427304797, 112.3783179235179, GalleryHeight, meshSizeCore};
Point(1168) = {-155.33432427304797, 102.62839277053718, GalleryHeight, meshSizeCore};
Point(1169) = {-194.93432427401422, -168.5718436865136, GalleryHeight, meshSizeCore};
Point(1170) = {-194.93432427401422, -163.171843685559, GalleryHeight, meshSizeCore};
Point(1171) = {-194.93432427401422, 36.82815631350968, GalleryHeight, meshSizeCore};
Point(1172) = {-194.93432427401422, 42.228156314464286, GalleryHeight, meshSizeCore};
Point(1173) = {-155.33432427304797, -168.57190000952687, GalleryHeight, meshSizeCore};
Point(1174) = {-155.33432427304797, -163.17189997050446, GalleryHeight, meshSizeCore};
Point(1175) = {-155.33432427304797, 36.828099991427734, GalleryHeight, meshSizeCore};
Point(1176) = {-155.33432427304797, 42.22809999145102, GalleryHeight, meshSizeCore};
Point(1178) = {-160.5343242760282, -168.57189261447638, GalleryHeight, meshSizeCore};
Point(1179) = {-160.5343242760282, -163.1718925795285, GalleryHeight, meshSizeCore};
Point(1180) = {-160.5343242760282, 36.82810738647822, GalleryHeight, meshSizeCore};
Point(1181) = {-160.5343242760282, 42.22810738743283, GalleryHeight, meshSizeCore};
Point(1185) = {-115.23432427900843, -233.97, GalleryHeight, meshSizeCore};
Point(1186) = {147.4847248259466, 62.22849326150026, GalleryHeight, meshSizeCloseToSeam};
Point(1187) = {197.4847140919882, 62.22849326150026, GalleryHeight, meshSizeCore};
Point(1188) = {147.48472712998046, -148.5720208064885, GalleryHeight, meshSizeCloseToSeam};
Point(1189) = {197.4847140919882, -148.5720208064885, GalleryHeight, meshSizeCore};


Point(1190) = {-115.23432427900843, 56.828797810478136, SeamHeight, meshSizeSeam};
Point(1191) = {-60.4151782640256, 72.22858866746537, SeamHeight, meshSizeSeam};
Point(1192) = {-44.35527517402079, 56.828797810478136, SeamHeight, meshSizeSeam};
Point(1193) = {147.4847248259466 + 4., 56.82851888751611, SeamHeight, meshSizeSeam};
Point(1194) = {197.48472482495708, 56.82837679947261, SeamHeight, meshSizeCore};
Point(1195) = {241.06384326395346, 56.82842257351149, SeamHeight, meshSizeCore};
Point(1196) = {241.06384326395346, -143.17125472647604, SeamHeight, meshSizeCore};
Point(1197) = {197.48471409198828, -143.17125238850713, SeamHeight, meshSizeCore};
Point(1198) = {147.48472712998046 + 4., -143.17124071449507, SeamHeight, meshSizeSeam};
Point(1199) = {-115.23432427900843, -143.17110019456595, SeamHeight, meshSizeSeam};

Point(1201) = {-230.0899094590568, 112.3783179235179, -BaseDepth, meshSizeRemoteCore};
Point(1202) = {241.06384326395346, 112.378317923517, -BaseDepth, meshSizeRemoteCore};
Point(1203) = {241.06384326395346, -233.97, -BaseDepth, meshSizeRemoteCore};
Point(1204) = {-230.0899094590568, -233.97, -BaseDepth, meshSizeRemoteCore};
Point(1205) = {-230.0899094590568, 112.3783179235179, CoreHeight, meshSizeRemoteCore};
Point(1206) = {241.06384326395346, 112.378317923517, CoreHeight, meshSizeRemoteCore};
Point(1207) = {241.06384326395346, -233.97, CoreHeight, meshSizeCore};
Point(1208) = {-230.0899094590568, -233.97, CoreHeight, meshSizeRemoteCore};


Line(143) = {226,821};
Line(144) = {312,907};
Line(145) = {314,909};
Line(146) = {318,913};
Line(147) = {319,914};
Line(148) = {320,915};
Line(149) = {321,916};
Line(150) = {400,995};
Line(151) = {401,996};
Line(152) = {403,998};
Line(153) = {404,999};
Line(154) = {407,1002};
Line(155) = {411,1006};
Line(156) = {415,1010};
Line(157) = {420,1015};
Line(158) = {421,1016};
Line(159) = {426,1021};
Line(160) = {427,1022};
Line(161) = {428,1023};
Line(162) = {432,1027};
Line(163) = {434,1029};
Line(164) = {436,1031};
Line(165) = {438,1033};
Line(166) = {440,1035};
Line(167) = {444,1039};
Line(168) = {445,1040};
Line(169) = {446,1041};
Line(170) = {450,1045};
Line(171) = {451,1046};
Line(172) = {455,1050};
Line(173) = {459,1054};
Line(174) = {460,1055};
Line(175) = {461,1056};
Line(176) = {462,1057};
Line(177) = {463,1058};
Line(178) = {464,1059};
Line(179) = {465,1060};
Line(180) = {466,1061};
Line(181) = {471,1066};
Line(182) = {472,1067};
Line(183) = {473,1068};
Line(184) = {474,1069};
Line(185) = {475,1070};
Line(186) = {476,1071};
Line(187) = {477,1072};
Line(188) = {478,1073};
Line(189) = {479,1074};
Line(190) = {484,1079};
Line(191) = {485,1080};
Line(192) = {486,1081};
Line(193) = {487,1082};
Line(194) = {488,1083};
Line(195) = {489,1084};
Line(196) = {490,1085};
Line(197) = {497,1092};
Line(198) = {498,1093};
Line(199) = {500,1095};
Line(200) = {501,1096};
Line(201) = {502,1097};
Line(202) = {504,1099};
Line(203) = {506,1101};
Line(204) = {507,1102};
Line(205) = {512,1107};
Line(206) = {513,1108};
Line(207) = {514,1109};
Line(208) = {515,1110};
Line(209) = {517,1112};
Line(210) = {518,1113};
Line(211) = {529,1124};
Line(212) = {530,1125};
Line(213) = {531,1126};
Line(214) = {532,1127};
Line(215) = {533,1128};
Line(216) = {535,1130};
Line(217) = {536,1131};
Line(218) = {537,1132};
Line(219) = {538,1133};
Line(220) = {539,1134};
Line(221) = {540,1135};
Line(222) = {541,1136};
Line(223) = {542,1137};
Line(224) = {543,1138};
Line(225) = {544,1139};
Line(226) = {545,1140};
Line(227) = {546,1141};
Line(228) = {548,1143};
Line(229) = {550,1145};
Line(230) = {559,1154};
Line(231) = {561,1156};
Line(232) = {562,1157};
Line(233) = {564,1159};
Line(234) = {565,1160};
Line(235) = {566,1161};
Line(236) = {567,1162};
Line(237) = {568,1163};
Line(238) = {569,1164};
Line(239) = {570,1165};
Line(240) = {571,1166};
Line(241) = {572,1167};
Line(242) = {573,1168};
Line(243) = {574,1169};
Line(244) = {575,1170};
Line(245) = {576,1171};
Line(246) = {577,1172};
Line(247) = {578,1173};
Line(248) = {579,1174};
Line(249) = {580,1175};
Line(250) = {581,1176};
Line(251) = {583,1178};
Line(252) = {584,1179};
Line(253) = {585,1180};
Line(254) = {586,1181};
Line(255) = {590,1185};
Line(256) = {591,1186};
Line(257) = {592,1187};
Line(258) = {593,1188};
Line(259) = {594,1189};
Line(260) = {1127,1128};
Line(261) = {1157,1107};
Line(262) = {1107,1163};
Line(263) = {1163,1162};
Line(264) = {1143,1160};
Line(265) = {1160,1161};
Line(266) = {1161,1156};
Line(267) = {1156,1143};
Line(268) = {1061,1159};
Line(269) = {1159,1181};
Line(270) = {1181,1172};
Line(271) = {1172,1061};
Line(272) = {1180,1171};
Line(273) = {1180,1139};
Line(274) = {1139,1140};
Line(275) = {1140,1171};
Line(276) = {1138,1131};
Line(277) = {1131,1130};
Line(278) = {1130,1056};
Line(279) = {1056,1138};
Line(280) = {1054,1133};
Line(281) = {1133,1179};
Line(282) = {1179,1170};
Line(283) = {1170,1054};
Line(284) = {1169,1178};
Line(285) = {1178,1132};
Line(286) = {1169,1154};
Line(287) = {1167,1168};
Line(288) = {1168,1080};
Line(289) = {1080,1079};
Line(290) = {1165,1067};
Line(291) = {1067,1066};
Line(292) = {1066,1166};
Line(293) = {1166,1165};
Line(294) = {1164,1069};
Line(295) = {1069,1071};
Line(296) = {1071,1176};
Line(297) = {1164,1176};
Line(298) = {1175,1070};
Line(299) = {1070,1068};
Line(300) = {1068,1141};
Line(301) = {1141,1175};
Line(302) = {1135,1082};
Line(303) = {1082,1081};
Line(304) = {1081,1134};
Line(305) = {1134,1135};
Line(306) = {1137,1023};
Line(307) = {1023,1085};
Line(308) = {1085,1174};
Line(309) = {1174,1137};
Line(310) = {1173,1084};
Line(311) = {1084,1083};
Line(312) = {1173,1136};
Line(313) = {1132,1136};
Line(314) = {1154,1128};
Line(315) = {1127,1157};
Line(316) = {1162,1167};
Line(317) = {1079,1073};
Line(318) = {1073,1072};
Line(319) = {1072,1095};
Line(320) = {1074,1093};
Line(321) = {1126,1112};
Line(322) = {1110,1145};
Line(323) = {1097,1096};
Line(324) = {1096,1095};
Line(325) = {1097,1101};
Line(326) = {1101,1102};
Line(327) = {1145,1109};
Line(328) = {1109,1108};
Line(329) = {1108,1110};
Line(330) = {1112,1113};
Line(331) = {1113,1125};
Line(332) = {1125,1126};
Line(333) = {1092,1093};
Line(334) = {1010,1092};
Line(335) = {1102,1099};
Line(336) = {1099,1124};
Line(337) = {1124,1058};
Line(338) = {821,1006};
Line(339) = {1006,1027};
Line(340) = {1010,1074};
Line(341) = {1027,1022};
Line(342) = {1060,1059};
Line(343) = {999,1057};
Line(344) = {1040,1039};
Line(345) = {1039,1035};
Line(346) = {1035,1015};
Line(347) = {1002,909};
Line(348) = {1050,1055};
Line(349) = {1046,1041};
Line(350) = {1022,1002};
Line(351) = {1021,1055};
Line(352) = {1029,1045};
Line(353) = {1021,1016};
Line(354) = {996,1033};
Line(355) = {995,1031};
Line(356) = {1031,1029};
Line(357) = {1033,1041};
Line(358) = {1046,996};
Line(359) = {1016,913};
Line(360) = {913,914};
Line(361) = {914,916};
Line(362) = {916,915};
Line(363) = {915,1050};
Line(364) = {821,907};
Line(365) = {998,1060};
Line(366) = {1057,1040};
Line(367) = {1083,1185};
Line(368) = {1015,1185};
Line(369) = {1058,1186};
Line(370) = {1186,1187};
Line(371) = {1187,1059};
Line(372) = {995,1188};
Line(373) = {1188,1189};
Line(374) = {1189,1045};
Line(375) = {1027, 1190};
Line(376) = {1190, 1199};
Line(377) = {1199, 1022};
Line(378) = {1199, 1198};
Line(379) = {1193, 1198};
Line(380) = {1198, 909};
Line(383) = {907, 1193};
Line(385) = {998, 1194};
Line(386) = {1194, 1193};
Line(387) = {1193, 1194};
Line(388) = {1194, 1197};
Line(389) = {1197, 999};
Line(391) = {1196, 1057};
Line(392) = {1196, 1197};
Line(393) = {1196, 1195};
Line(394) = {1194, 1195};
Line(395) = {1195, 1060};
Line(396) = {1193, 1192};
Line(397) = {1192, 821};
Line(398) = {1006, 1191};
Line(399) = {1191, 1192};
Line(400) = {1191, 1190};
Line(401) = {467, 1204};
Line(402) = {1204, 1201};
Line(403) = {1201, 468};
Line(404) = {468, 1205};
Line(405) = {1205, 1208};
Line(406) = {1208, 467};
Line(407) = {1204, 1203};
Line(408) = {1203, 589};
Line(409) = {589, 1207};
Line(410) = {1207, 1208};
Line(411) = {1207, 1206};
Line(412) = {1206, 588};
Line(413) = {588, 1202};
Line(414) = {1202, 1203};
Line(415) = {1206, 1205};
Line(416) = {1202, 1201};
//+
Curve Loop(31) = {8, 228, -267, -231};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {261, -205, -2, 232};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {262, -237, -3, 205};
//+
Plane Surface(33) = {33};
//+
Curve Loop(34) = {4, 236, -263, -237};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {29, 242, -287, -241};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {242, 288, -191, -30};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {187, -318, -188, 67};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {289, -190, -31, 191};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {319, -199, -69, 187};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {200, -323, -201, 73};
//+
Plane Surface(40) = {40};
//+
Curve Loop(41) = {325, -203, -75, 201};
//+
Plane Surface(41) = {41};
//+
Curve Loop(42) = {326, -204, -76, 203};
//+
Plane Surface(42) = {42};
//+
Curve Loop(43) = {336, -211, -88, 202};
//+
Plane Surface(43) = {43};
//+
Curve Loop(44) = {327, -207, -77, 229};
//+
Plane Surface(44) = {44};
//+
Curve Loop(45) = {322, -229, -72, 208};
//+
Plane Surface(45) = {45};
//+
Curve Loop(46) = {328, -206, -78, 207};
//+
Plane Surface(46) = {46};
//+
Curve Loop(47) = {329, -208, -79, 206};
//+
Plane Surface(47) = {47};
//+
Curve Loop(48) = {71, 209, -321, -213};
//+
Plane Surface(48) = {48};
//+
Curve Loop(49) = {80, 210, -330, -209};
//+
Plane Surface(49) = {49};
//+
Curve Loop(50) = {81, 212, -331, -210};
//+
Plane Surface(50) = {50};
//+
Curve Loop(51) = {332, -213, -82, 212};
//+
Plane Surface(51) = {51};
//+
Curve Loop(52) = {83, 198, -333, -197};
//+
Plane Surface(52) = {52};
//+
Curve Loop(53) = {70, 198, -320, -189};
//+
Plane Surface(53) = {53};
//+
Curve Loop(54) = {95, 189, -340, -156};
//+
Curve Loop(55) = {337, -177, -89, 211};
//+
Plane Surface(54) = {55};
//+
Curve Loop(56) = {137, 256, -369, -177};
//+
Plane Surface(55) = {56};
//+
Curve Loop(57) = {338, 398, 399, 397};
//+
Plane Surface(56) = {57};
//+
Curve Loop(58) = {397, 364, 383, 396};
//+
Plane Surface(57) = {58};
//+
Curve Loop(59) = {379, 380, -145, -100, 144, 383};
//+
Plane Surface(58) = {59};
//+
Line(417) = {998, 907};
//+
Curve Loop(60) = {378, -379, 396, -399, 400, 376};
//+
Plane Surface(59) = {60};
//+
Curve Loop(61) = {84, 197, -334, -156};
//+
Plane Surface(60) = {61};
//+
Plane Surface(61) = {54};
//+
Curve Loop(62) = {339, 375, -400, -398};
//+
Plane Surface(62) = {62};
//+
Curve Loop(63) = {94, 162, -339, -155};
//+
Plane Surface(63) = {63};
//+
Curve Loop(64) = {93, 155, -338, -143};
//+
Plane Surface(64) = {64};
//+
Curve Loop(65) = {128, 144, -364, -143};
//+
Plane Surface(65) = {65};
//+
Curve Loop(66) = {102, 153, -389, -388, -385, -152};
//+
Plane Surface(66) = {66};
//+
Curve Loop(67) = {343, -391, 392, 389};
//+
Plane Surface(67) = {67};
//+
Curve Loop(68) = {138, 257, -370, -256};
//+
Plane Surface(68) = {68};
//+
Curve Loop(69) = {139, 178, -371, -257};
//+
Plane Surface(69) = {69};
//+
Curve Loop(70) = {417, 383, -386, -385};
//+
Plane Surface(70) = {70};
//+
Curve Loop(71) = {365, -395, -394, -385};
//+
Plane Surface(71) = {71};
//+
Curve Loop(72) = {291, -181, -33, 182};
//+
Plane Surface(72) = {72};
//+
Curve Loop(73) = {290, -182, -32, 239};
//+
Plane Surface(73) = {73};
//+
Curve Loop(74) = {293, -239, -35, 240};
//+
Plane Surface(74) = {74};
//+
Curve Loop(75) = {34, 240, -292, -181};
//+
Plane Surface(75) = {75};
//+
Curve Loop(76) = {264, -234, -5, 228};
//+
Plane Surface(76) = {76};
//+
Curve Loop(77) = {6, 235, -265, -234};
//+
Plane Surface(77) = {77};
//+
Curve Loop(78) = {294, -184, -36, 238};
//+
Line(418) = {999, 909};
//+
Line(419) = {1197, 1198};
//+
Curve Loop(79) = {388, 419, -379, -386};
//+
Plane Surface(78) = {79};
//+
Curve Loop(80) = {101, 176, -391, 393, 395, -179};
//+
Plane Surface(79) = {80};

Curve Loop(81) = {378, 380, -347, -350, -377};
//+
Plane Surface(80) = {81};
//+
Curve Loop(82) = {160, 350, -154, -112};
//+
Plane Surface(81) = {82};
//+
Curve Loop(83) = {109, 145, -347, -154};
//+
Plane Surface(82) = {83};
//+
Curve Loop(84) = {418, -380, -419, 389};
//+
Plane Surface(83) = {84};
//+
Curve Loop(85) = {140, 258, -372, -150};
//+
Plane Surface(84) = {85};
//+
Curve Loop(86) = {141, 259, -373, -258};
//+
Plane Surface(85) = {86};
//+
Curve Loop(87) = {259, 374, -170, -142};
//+
Plane Surface(86) = {87};
//+
Curve Loop(88) = {352, -170, -115, 163};
//+
Plane Surface(87) = {88};
//+
Curve Loop(89) = {343, -176, -103, 153};
//+
Plane Surface(88) = {89};
//+
Curve Loop(90) = {104, 167, -344, -168};
//+
Plane Surface(89) = {90};
//+
Curve Loop(91) = {345, -166, -105, 167};
//+
Plane Surface(90) = {91};
//+
Curve Loop(92) = {106, 157, -346, -166};
//+
Plane Surface(91) = {92};
//+
Curve Loop(93) = {120, 163, -356, -164};
//+
Plane Surface(92) = {93};
//+
Curve Loop(94) = {355, -164, -118, 150};
//+
Plane Surface(93) = {94};
//+
Curve Loop(95) = {165, 357, -169, -121};
//+
Plane Surface(94) = {95};
//+
Curve Loop(96) = {117, 165, -354, -151};
//+
Plane Surface(95) = {96};
//+
Curve Loop(97) = {122, 151, -358, -171};
//+
Plane Surface(96) = {97};
//+
Curve Loop(98) = {349, -169, -111, 171};
//+
Plane Surface(97) = {98};
//+
Curve Loop(99) = {348, -174, -110, 172};
//+
Plane Surface(98) = {99};
//+
Curve Loop(100) = {113, 174, -351, -159};
//+
Plane Surface(99) = {100};
//+
Curve Loop(101) = {353, -158, -116, 159};
//+
Plane Surface(100) = {101};
//+
Curve Loop(102) = {123, 146, -359, -158};
//+
Plane Surface(101) = {102};
//+
Curve Loop(103) = {124, 147, -360, -146};
//+
Plane Surface(102) = {103};
//+
Curve Loop(104) = {126, 148, -362, -149};
//+
Plane Surface(103) = {104};
//+
Curve Loop(105) = {361, -149, -125, 147};
//+
Plane Surface(104) = {105};
//+
Curve Loop(106) = {127, 172, -363, -148};
//+
Plane Surface(105) = {106};
//+
Curve Loop(107) = {134, 255, -368, -157};
//+
Plane Surface(106) = {107};
//+
Curve Loop(108) = {307, -196, -49, 161};
//+
Plane Surface(107) = {108};
//+
Curve Loop(109) = {306, -161, -48, 223};
//+
Plane Surface(108) = {109};
//+
Curve Loop(110) = {308, -248, -50, 196};
//+
Plane Surface(109) = {110};
//+
Curve Loop(111) = {309, -223, -51, 248};
//+
Curve Loop(112) = {377, -341, 375, 376};
//+
Plane Surface(110) = {112};
//+
Curve Loop(113) = {160, -341, -162, 96};
//+
Plane Surface(111) = {113};
//+
Curve Loop(114) = {45, 192, -303, -193};
//+
Plane Surface(112) = {114};
//+
Curve Loop(115) = {46, 220, -304, -192};
//+
Plane Surface(113) = {115};
//+
Curve Loop(116) = {277, -216, -18, 217};
//+
Plane Surface(114) = {116};
//+
Curve Loop(117) = {44, 193, -302, -221};
//+
Plane Surface(115) = {117};
//+
Curve Loop(118) = {305, -221, -47, 220};
//+
Plane Surface(116) = {118};
//+
Curve Loop(119) = {19, 175, -278, -216};
//+
Plane Surface(117) = {119};
//+
Curve Loop(120) = {20, 224, -279, -175};
//+
Plane Surface(118) = {120};
//+
Curve Loop(121) = {17, 217, -276, -224};
//+
Plane Surface(119) = {121};
//+
Curve Loop(122) = {1, 215, -260, -214};
//+
Plane Surface(120) = {122};
//+
Curve Loop(123) = {11, 246, -270, -254};
//+
Plane Surface(121) = {123};
//+
Curve Loop(124) = {10, 254, -269, -233};
//+
Plane Surface(122) = {124};
//+
Curve Loop(125) = {9, 233, -268, -180};
//+
Plane Surface(123) = {125};
//+
Curve Loop(126) = {12, 180, -271, -246};
//+
Plane Surface(124) = {126};
//+
Curve Loop(127) = {250, -297, -238, 39};
//+
Curve Loop(128) = {37, 186, -295, -184};
//+
Plane Surface(125) = {128};
//+
Curve Loop(129) = {38, 250, -296, -186};
//+
Plane Surface(126) = {129};
//+
Plane Surface(127) = {127};
//+
Curve Loop(130) = {13, 245, -272, -253};
//+
Plane Surface(128) = {130};
//+
Curve Loop(131) = {16, 245, -275, -226};
//+
Plane Surface(129) = {131};
//+
Curve Loop(132) = {274, -226, -15, 225};
//+
Plane Surface(130) = {132};
//+
Curve Loop(133) = {14, 225, -273, -253};
//+
Plane Surface(131) = {133};
//+
Curve Loop(134) = {43, 249, -301, -227};
//+
Plane Surface(132) = {134};
//+
Curve Loop(135) = {40, 185, -298, -249};
//+
Plane Surface(133) = {135};
//+
Curve Loop(136) = {41, 183, -299, -185};
//+
Plane Surface(134) = {136};
//+
Curve Loop(137) = {300, -227, -42, 183};
//+
Plane Surface(135) = {137};
//+
Curve Loop(138) = {284, -251, -26, 243};
//+
Plane Surface(136) = {138};
//+
Curve Loop(139) = {244, 283, -173, -25};
//+
Curve Loop(140) = {21, 219, -280, -173};
//+
Plane Surface(137) = {140};
//+
Plane Surface(138) = {139};
//+
Curve Loop(141) = {282, -244, -23, 252};
//+
Plane Surface(139) = {141};
//+
Curve Loop(142) = {22, 252, -281, -219};
//+
Plane Surface(140) = {142};
//+
Curve Loop(143) = {28, 230, -286, -243};
//+
Plane Surface(141) = {143};
//+
Curve Loop(144) = {27, 218, -285, -251};
//+
Plane Surface(142) = {144};
//+
Curve Loop(145) = {55, 222, -312, -247};
//+
Plane Surface(143) = {145};
//+
Curve Loop(146) = {52, 195, -310, -247};
//+
Plane Surface(144) = {146};
//+
Curve Loop(147) = {53, 194, -311, -195};
//+
Plane Surface(145) = {147};
//+
Curve Loop(148) = {394, -393, 392, -388};
//+
Plane Surface(146) = {148};
//+
Curve Loop(149) = {260, -314, -286, 284, 285, 313, -312, 310, 311, 367, -368, -346, -345, -344, -366, -343, 418, -347, -350, -341, -339, -338, 364, -417, 365, 342, -371, -370, -369, -337, -336, -335, -326, -325, 323, 324, -319, -318, -317, -289, -288, -287, -316, -263, -262, -261, -315};
//+
Curve Loop(150) = {264, 265, 266, 267};
//+
Curve Loop(151) = {292, 293, 290, 291};
//+
Curve Loop(152) = {334, 333, -320, -340};
//+
Curve Loop(153) = {332, 321, 330, 331};
//+
Curve Loop(154) = {329, 322, 327, 328};
//+
Curve Loop(155) = {268, 269, 270, 271};
//+
Curve Loop(156) = {296, -297, 294, 295};
//+
Curve Loop(157) = {272, -275, -274, -273};
//+
Curve Loop(158) = {300, 301, 298, 299};
//+
Curve Loop(159) = {276, 277, 278, 279};
//+
Curve Loop(160) = {305, 302, 303, 304};
//+
Curve Loop(161) = {280, 281, 282, 283};
//+
Curve Loop(162) = {308, 309, 306, 307};
//+
Curve Loop(163) = {351, -348, -363, -362, -361, -360, -359, -353};
//+
Curve Loop(164) = {349, -357, -354, -358};
//+
Curve Loop(165) = {355, 356, 352, -374, -373, -372};
//+
Plane Surface(147) = {149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165};
//+
Curve Loop(166) = {406, -59, 215, -314, -230, -57, 218, 313, -222, -54, 194, 367, -255, 132, 409, 410};
//+
Plane Surface(148) = {166};
//+
Curve Loop(167) = {133, 168, -366, -391, 393, 395, 342, -178, 91, 202, -335, -204, 92, 200, 324, -199, 136, -412, -411, -409};
//+
Plane Surface(149) = {167};
//+
Curve Loop(168) = {408, -132, -131, 54, -56, 57, 58, 59, 401, 407};
//+
Plane Surface(150) = {168};
//+
Curve Loop(169) = {133, -130, -101, 99, 91, -87, 92, 74, 136, 413, 414, 408};
//+
Plane Surface(151) = {169};
//+
Curve Loop(170) = {406, 60, 404, 405};
//+
Plane Surface(152) = {170};
//+
Curve Loop(171) = {60, -403, -402, -401};
//+
Plane Surface(153) = {171};
//+
Curve Loop(172) = {405, -410, 411, 415};
//+
Curve Loop(173) = {414, -407, 402, -416};
//+
Plane Surface(154) = {172};
//+
Plane Surface(155) = {173};
//+
Curve Loop(174) = {403, 61, 62, 63, 64, 65, 66, 135, 413, 416};
//+
Plane Surface(156) = {174};
//+
Curve Loop(175) = {415, -404, 61, 214, 315, -232, 63, 236, 316, -241, 65, 190, 317, -188, 135, -412};
//+
Plane Surface(157) = {175};
//+
Delete {
  Surface{147}; 
}
//+
Delete {
  Surface{30}; 
}
//+
Delete {
  Surface{80}; Surface{81}; Surface{82}; Surface{21}; 
}
//+
Delete {
  Curve{350}; Curve{154}; Curve{347}; Curve{109}; Curve{112}; 
}
//+
Delete {
  Point{1002}; Point{407}; 
}
//+
Line(420) = {1022, 909};
//+
Line(421) = {427, 314};
//+
Curve Loop(176) = {160, 420, -145, -421};
//+
Plane Surface(158) = {176};
//+
Curve Loop(177) = {420, -380, -378, 377};
//+
Plane Surface(159) = {177};
//+
Delete {
  Surface{110}; Surface{111}; 
}
//+
Delete {
  Surface{63}; Surface{62}; 
}
//+
Delete {
  Surface{64}; Surface{56}; 
}
//+
Delete {
  Surface{65}; 
}
//+
Delete {
  Surface{57}; 
}
//+
Delete {
  Curve{341}; Curve{96}; 
}
//+
Delete {
  Curve{128}; 
}
//+
Delete {
  Curve{364}; 
}
//+
Delete {
  Curve{93}; 
}
Point(1210) = {-49.9861, 62.2293, 0, meshSizeCloseToSeam};
Point(1212) = {-49.9861, 62.2293, GalleryHeight, meshSizeCloseToSeam};

Point(1209) = {-115.2343242790084, 62.226, 0, meshSizeCloseToSeam};
Point(1211) = {-115.2343242790084, 62.2293, GalleryHeight, meshSizeCloseToSeam};
Point(1214) = {-115.23432427900843, 56.82856875646394, 0., meshSizeSeam};
Point(1215) = {-115.23432427900843, 56.82856875646394, GalleryHeight, meshSizeSeam};

//Point(1069) = {-120.63432427705266, 56.82856875646394, GalleryHeight, meshSizeCloseToSeam};//+
Recursive Delete {
  Curve{396}; Curve{399}; Curve{397}; Curve{338}; 
}
//+
Recursive Delete {
  Curve{396}; 
}
//+
Recursive Delete {
  Curve{396}; Curve{399}; 
}
//+
Delete {
  Surface{59}; 
}
//+
Line(422) = {1215, 1022};
//+
Line(423) = {1214, 1215};
//+
Line(424) = {1214, 427};
//+
Line(425) = {1209, 1211};
//+
Line(426) = {1027, 1211};
//+
Line(427) = {432, 1209};
//+
Line(428) = {1211, 1215};
//+
Curve Loop(178) = {427, 425, -426, -162};
//+
Plane Surface(160) = {178};
//+
Curve Loop(179) = {428, 422, -377, -376, -375, 426};
//+
Plane Surface(161) = {179};
//+
Curve Loop(180) = {424, 160, -422, -423};
//+
Plane Surface(162) = {180};
//+
Plane Surface(163) = {63};
//+
Plane Surface(164) = {62};
//+
Delete {
  Surface{164}; 
}
//+
Delete {
  Surface{161}; 
}
//+
Delete {
  Curve{398}; Curve{400}; Curve{375}; 
}
//+
Delete {
  Curve{422}; 
}
//+
Line(429) = {1209, 1210};
//+
Line(430) = {1210, 1212};
//+
Line(431) = {1210, 411};
//+
Line(432) = {1212, 1006};
//+
Line(433) = {1212, 1211};
//+
Delete {
  Curve{396}; Curve{399}; Curve{143}; 
}
//+
Line(434) = {1193, 1190};
//+
Line(435) = {1190, 1215};
//+
Line(436) = {1215, 907};
//+
Line(437) = {1214, 312};
//+
Delete {
  Curve{428}; 
}
//+
Delete {
  Point{1191}; Point{1192}; Point{821}; Point{226}; 
}
//+
Curve Loop(181) = {431, 155, -432, -430};
Plane Surface(164) = {181};

Curve Loop(182) = {429, 430, 433, -425};
Plane Surface(165) = {182};

Curve Loop(183) = {437, 144, -436, -423};
Plane Surface(166) = {183};

Curve Loop(184) = {434, 435, 436, 383};
Plane Surface(167) = {184};

Curve Loop(185) = {422, -377, -376, 435};
Plane Surface(168) = {185};

Curve Loop(186) = {429, 431, 94, 427};
Plane Surface(169) = {186};
Plane Surface(170) = {78};

Curve Loop(187) = {424, 421, -100, -437};
Plane Surface(171) = {187};

Curve Loop(188) = {376, 378, -379, 434};
Plane Surface(172) = {188};

//+
Recursive Delete {
  Curve{359}; Curve{123}; Curve{146}; Curve{360}; Curve{124}; Curve{147}; Curve{125}; Curve{361}; Curve{149}; Curve{362}; Curve{126}; Curve{148}; Curve{363}; Curve{127}; 
}
//+
Recursive Delete {
  Surface{27}; 
}
//+
Delete {
  Surface{101}; Surface{102}; Surface{104}; Surface{103}; Surface{105}; 
}
//+
Delete {
  Curve{123}; Curve{359}; Curve{146}; Curve{360}; Curve{124}; Curve{147}; Curve{361}; Curve{125}; Curve{362}; Curve{149}; Curve{126}; Curve{148}; Curve{363}; Curve{127}; 
}
//+
Delete {
  Point{914}; Point{916}; Point{319}; Point{321}; Point{913}; Point{915}; Point{318}; Point{320}; 
}
//+
Line(438) = {1016, 1050};
//+
Line(439) = {421, 455};
//+
Curve Loop(189) = {439, 172, -438, -158};
//+
Plane Surface(173) = {189};
//+
Curve Loop(190) = {439, 110, -113, 116};
//+
Plane Surface(174) = {190};
//+
Curve Loop(191) = {315, 261, 262, 263, 316, 287, 288, 289, 317, 318, 319, -324, -323, 325, 326, 335, 336, 337, 369, 370, 371, -342, -365, 417, -436, 422, 420, -418, 343, 366, 344, 345, 346, 368, -367, -311, -310, 312, -313, -285, -284, 286, 314, -260};
//+
Curve Loop(192) = {339, 426, -433, 432};
//+
Curve Loop(193) = {351, -348, -438, -353};
//+
Plane Surface(175) = {150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 164, 165, 191, 192, 193};
//+
Plane Surface(176) = {3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 17, 18, 19, 26, 29, 30, 186, 187, 190};
Point(1317) = {-65.5792996730306, -178.422046875465, 1.5, meshSizeElectrodes};
Point{1317} In Surface {173};
Point(1318) = {-60.5803996730247, -178.422046875465, 1.5, meshSizeElectrodes};
Point{1318} In Surface {173};
Point(1319) = {-52.3857996730367, -175.513446875499, 1.5, meshSizeElectrodes};
Point{1319} In Surface {98};
Point(1320) = {-48.7492996730143, -171.876946875476, 1.5, meshSizeElectrodes};
Point{1320} In Surface {98};
Point(1321) = {-45.1127996730502, -168.24044687557, 1.5, meshSizeElectrodes};
Point{1321} In Surface {98};
Point(1322) = {-41.4761996730231, -164.603946875548, 1.5, meshSizeElectrodes};
Point{1322} In Surface {98};
Point(1323) = {-37.8396996730589, -160.967446875526, 1.5, meshSizeElectrodes};
Point{1323} In Surface {98};
Point(1324) = {-34.2031996730366, -157.330946875503, 1.5, meshSizeElectrodes};
Point{1324} In Surface {98};
Point(1325) = {-30.5666996730142, -153.694446875481, 1.5, meshSizeElectrodes};
Point{1325} In Surface {98};
Point(1326) = {-26.9301996730501, -150.057946875459, 1.5, meshSizeElectrodes};
Point{1326} In Surface {98};
Point(1327) = {-21.5133996730438, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1327} In Surface {158};
Point(1328) = {-16.3940996730234, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1328} In Surface {158};
Point(1329) = {-11.2727996730246, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1329} In Surface {158};
Point(1330) = {-6.15159967303043, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1330} In Surface {158};
Point(1331) = {-1.0303996730363, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1331} In Surface {158};
Point(1332) = {4.09080032695783, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1332} In Surface {158};
Point(1333) = {9.21200032695197, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1333} In Surface {158};
Point(1334) = {14.3332003269461, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1334} In Surface {158};
Point(1335) = {19.4544003269402, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1335} In Surface {158};
Point(1336) = {24.5756003269926, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1336} In Surface {158};
Point(1337) = {29.6968003269867, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1337} In Surface {158};
Point(1338) = {34.8181003269856, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1338} In Surface {158};
Point(1339) = {39.9393003269797, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1339} In Surface {158};
Point(1340) = {45.0605003269739, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1340} In Surface {158};
Point(1341) = {50.181700326968, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1341} In Surface {158};
Point(1342) = {55.3029003269621, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1342} In Surface {158};
Point(1343) = {60.4241003269563, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1343} In Surface {158};
Point(1344) = {65.5453003269504, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1344} In Surface {158};
Point(1345) = {70.6665003269445, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1345} In Surface {158};
Point(1346) = {75.7878003269434, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1346} In Surface {158};
Point(1347) = {80.9090003269375, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1347} In Surface {158};
Point(1348) = {86.0302003269899, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1348} In Surface {158};
Point(1349) = {91.151400326984, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1349} In Surface {158};
Point(1350) = {96.2726003269781, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1350} In Surface {158};
Point(1351) = {101.393800326972, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1351} In Surface {158};
Point(1352) = {106.515000326966, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1352} In Surface {158};
Point(1353) = {111.636200326961, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1353} In Surface {158};
Point(1354) = {116.757500326959, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1354} In Surface {158};
Point(1355) = {121.878700326954, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1355} In Surface {158};
Point(1356) = {126.999900326948, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1356} In Surface {158};
Point(1357) = {132.121100326942, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1357} In Surface {158};
Point(1358) = {137.242300326994, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1358} In Surface {158};
Point(1359) = {142.363500326988, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1359} In Surface {158};
Point(1360) = {147.484700326982, -143.171246875543, 1.5, meshSizeElectrodes};
Point{1360} In Surface {158};
Point(1423) = {-45.0155996730318, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1423} In Surface {166};
Point(1424) = {-40.3192996730213, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1424} In Surface {166};
Point(1425) = {-35.6241996730096, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1425} In Surface {166};
Point(1426) = {-30.929099673056, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1426} In Surface {166};
Point(1427) = {-26.2339996730443, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1427} In Surface {166};
Point(1428) = {-21.5388996730326, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1428} In Surface {166};
Point(1429) = {-16.8437996730208, 56.8283531244379, 1.5, meshSizeElectrodes};
Point{1429} In Surface {166};
Point(1430) = {-12.1486996730091, 56.8288531245198, 1.5, meshSizeElectrodes};
Point{1430} In Surface {166};
Point(1431) = {-7.45359967305558, 56.8288531245198, 1.5, meshSizeElectrodes};
Point{1431} In Surface {166};
Point(1432) = {-2.75849967304384, 56.8288531245198, 1.5, meshSizeElectrodes};
Point{1432} In Surface {166};
Point(1433) = {1.93660032696789, 56.8288531245198, 1.5, meshSizeElectrodes};
Point{1433} In Surface {166};
Point(1434) = {6.63170032697963, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1434} In Surface {166};
Point(1435) = {11.3268003269914, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1435} In Surface {166};
Point(1436) = {16.0219003269449, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1436} In Surface {166};
Point(1437) = {20.7170003269566, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1437} In Surface {166};
Point(1438) = {25.4121003269684, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1438} In Surface {166};
Point(1439) = {30.1072003269801, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1439} In Surface {166};
Point(1440) = {34.8023003269918, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1440} In Surface {166};
Point(1441) = {39.4974003269454, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1441} In Surface {166};
Point(1442) = {44.1925003269571, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1442} In Surface {166};
Point(1443) = {48.8876003269688, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1443} In Surface {166};
Point(1444) = {53.5827003269806, 56.8288531245198, 1.5, meshSizeElectrodes};
Point{1444} In Surface {166};
Point(1445) = {58.2778003269923, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1445} In Surface {166};
Point(1446) = {62.9729003269458, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1446} In Surface {166};
Point(1447) = {67.6680003269576, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1447} In Surface {166};
Point(1448) = {72.3631003269693, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1448} In Surface {166};
Point(1449) = {77.058200326981, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1449} In Surface {166};
Point(1450) = {81.7533003269928, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1450} In Surface {166};
Point(1451) = {86.4484003269463, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1451} In Surface {166};
Point(1452) = {91.143500326958, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1452} In Surface {166};
Point(1453) = {95.8386003269698, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1453} In Surface {166};
Point(1454) = {100.533700326981, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1454} In Surface {166};
Point(1455) = {105.228800326993, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1455} In Surface {166};
Point(1456) = {109.923900326947, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1456} In Surface {166};
Point(1457) = {114.619000326958, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1457} In Surface {166};
Point(1458) = {119.31410032697, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1458} In Surface {166};
Point(1459) = {124.009200326982, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1459} In Surface {166};
Point(1460) = {128.704300326994, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1460} In Surface {166};
Point(1461) = {133.399400326947, 56.8286531245103, 1.5, meshSizeElectrodes};
Point{1461} In Surface {166};
Point(1462) = {138.094500326959, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1462} In Surface {166};
Point(1463) = {142.789600326971, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1463} In Surface {166};
Point(1464) = {147.484700326982, 56.8285531244474, 1.5, meshSizeElectrodes};
Point{1464} In Surface {166};
//+
Surface Loop(1) = {171, 172, 168, 159, 167, 58, 158, 162, 166};
//+
Volume(1) = {1};
//+
Curve Loop(194) = {129, 179, -365, -152};
//+
Plane Surface(177) = {194};
//+
Surface Loop(2) = {23, 146, 71, 67, 79, 88, 66, 177};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {155, 151, 176, 3, 10, 156, 153, 150, 1, 8, 9, 2, 15, 16, 17, 18, 28, 169, 14, 13, 12, 5, 11, 6, 7, 174, 24, 26, 29, 23, 20, 19, 171, 4};
//+
Volume(3) = {3};
//+
Curve Loop(195) = {266, -231, -7, 235};
//+
Plane Surface(178) = {195};
//+
Plane Surface(179) = {111};
//+
Surface Loop(4) = {154, 152, 148, 120, 175, 76, 77, 178, 31, 72, 75, 74, 73, 61, 53, 52, 60, 50, 51, 48, 49, 46, 47, 45, 44, 124, 123, 122, 121, 125, 126, 127, 170, 131, 130, 129, 128, 134, 135, 132, 133, 118, 119, 114, 117, 113, 116, 115, 112, 138, 139, 140, 137, 107, 109, 179, 108, 96, 95, 94, 97, 84, 85, 86, 87, 92, 93, 141, 136, 142, 143, 144, 145, 106, 91, 90, 89, 149, 69, 68, 55, 54, 43, 42, 41, 40, 39, 37, 157, 32, 33, 34, 35, 36, 38, 83, 78, 70, 164, 163, 160, 165, 100, 173, 98, 99, 1, 2, 3, 4, 5, 6, 12, 172, 168, 159, 167, 146, 71, 67, 17, 18, 19, 20, 28, 169, 16, 15, 14, 13, 11, 10, 9, 7, 8, 174, 26, 29, 24};
//+
Volume(4) = {4};

Point(1501) = {-230.0899094590568-PaddingX, 112.3783179235179+PaddingY, -BaseDepth-PaddingZ, meshSizePadding };
Point(1502) = {241.06384326395346+PaddingX, 112.378317923517+PaddingY, -BaseDepth-PaddingZ, meshSizePadding };
Point(1503) = {241.06384326395346+PaddingX, -233.97-PaddingY, -BaseDepth-PaddingZ, meshSizePadding };
Point(1504) = {-230.0899094590568-PaddingX, -233.97-PaddingY, -BaseDepth-PaddingZ, meshSizePadding };
Point(1505) = {-230.0899094590568-PaddingX, 112.3783179235179+PaddingY, CoreHeight+PaddingZ, meshSizePadding };
Point(1506) = {241.06384326395346+PaddingX, 112.378317923517+PaddingY, CoreHeight+PaddingZ, meshSizePadding };
Point(1507) = {241.06384326395346+PaddingX, -233.97-PaddingY, CoreHeight+PaddingZ, meshSizePadding };
Point(1508) = {-230.0899094590568-PaddingX, -233.97-PaddingY, CoreHeight+PaddingZ, meshSizePadding };//+
Line(440) = {1505, 1508};
//+
Line(441) = {1508, 1507};
//+
Line(442) = {1507, 1503};
//+
Line(443) = {1503, 1504};
//+
Line(444) = {1504, 1501};
//+
Line(445) = {1501, 1505};
//+
Line(446) = {1505, 1506};
//+
Line(447) = {1501, 1502};
//+
Line(448) = {1502, 1506};
//+
Line(449) = {1506, 1507};
//+
Line(450) = {1503, 1502};
//+
Line(451) = {1504, 1508};
//+
Curve Loop(196) = {58, 215, -314, -230};
//+
Plane Surface(180) = {196};
//+
Curve Loop(197) = {56, 222, -313, -218};
//+
Plane Surface(181) = {197};
//+
Curve Loop(198) = {367, -255, -131, 194};
//+
Plane Surface(182) = {198};
//+
Curve Loop(199) = {62, 232, -315, -214};
//+
Plane Surface(183) = {199};
//+
Curve Loop(200) = {64, 241, -316, -236};
//+
Plane Surface(184) = {200};
//+
Curve Loop(201) = {317, -188, -66, 190};
//+
Plane Surface(185) = {201};
//+
Curve Loop(202) = {366, -168, -130, 176};
//+
Plane Surface(186) = {202};
//+
Curve Loop(203) = {342, -178, -99, 179};
//+
Plane Surface(187) = {203};
//+
Curve Loop(204) = {335, -202, -87, 204};
//+
Plane Surface(188) = {204};
//+
Curve Loop(205) = {74, 199, -324, -200};
//+
Plane Surface(189) = {205};
//+
Curve Loop(206) = {449, 442, 450, 448};
//+
Plane Surface(190) = {206};
//+
Curve Loop(207) = {441, 442, 443, 451};
//+
Plane Surface(191) = {207};
//+
Curve Loop(208) = {444, 445, 440, -451};
//+
Plane Surface(192) = {208};
//+
Curve Loop(209) = {449, -441, -440, 446};
//+
Plane Surface(193) = {209};
//+
Curve Loop(210) = {447, 448, -446, -445};
//+
Plane Surface(194) = {210};
//+
Curve Loop(211) = {450, -447, -444, -443};
//+
Plane Surface(195) = {211};

//+
Surface Loop(5) = {191, 193, 190, 195, 194, 192};
//+
Surface Loop(6) = {153, 156, 151, 155, 150, 154, 152, 148, 79, 186, 188, 187, 189, 183, 184, 185, 157, 149, 180, 181, 182};
//+
Volume(5) = {5, 6};
//+
Physical Surface("Faces", 452) = {192, 193, 194, 190, 191, 195};
//+
Physical Volume("Goaf", 453) = {2};
Physical Volume("Seam", 454) = {1};
Physical Volume("Core", 455) = {4};
Physical Volume("Base", 456) = {3};
Physical Volume("Padding", 457) = {5};

