#!/usr/bin/python3
import numpy as np
import importlib, sys, os
import argparse

parser = argparse.ArgumentParser(description='converts a ZZ ERT file and creates fingal inp[ut files.', epilog="l.gross@uq.edu.au, version 12/01/2019")
parser.add_argument(dest='zzfile', metavar='zzfile', type=str, help='ZZ data file')
parser.add_argument('--NaN', '-N',  dest='NaN', default="9", type=str, help="string indicating faulty data")
parser.add_argument('--datafile' '-D',  dest='filename', type=str, help="file name of data file extracked. if not set survey name + extension 'csv' is used.")
parser.add_argument('--wenner', '-w',  dest='wenner', action='store_true', default=False, help="generates a Wenner survey otherwise a full dipole-dipole survey is created.")
parser.add_argument('--nogeo', '-n',  dest='nogeo', action='store_true', default=False, help="if set no geo file for gmsh is generated")
parser.add_argument('--geofile', '-G',  dest='geofile', type=str,  help="name of the geofile. If not set the survey name + extension 'geo' is used.")
parser.add_argument('--stationfile', '-S',  dest='stationfile', type=str,  help="name of the station/electrode file in format (A,B,M,N). If not set the survey name + extension 'loc' is used.")
parser.add_argument('--padding', '-p',  dest='relPadding', type=float, default=30, help="padding in %% relative to survey extend (forming the bounding box)")
parser.add_argument('--range', '-r',  dest='relRange', type=float, default=30, help="range of survey in %% relative to survey extend. The will be added to the extend and be updated during inversion")
parser.add_argument('--configfile', '-c',  dest='configfile', type=str, help="name of the configuarion file to be generated. If not set the survey name + extension 'Ypadding' is used.")
parser.add_argument('--centralize', '-z',  dest='centralize', action='store_true', default=False, help="the center of the survey is moved to (0,0).")

parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=False, help="shows more information.")
args = parser.parse_args()

f=open(args.zzfile ,'r')
ll=f.readline()
SurveyName=ll[ll.index(":",)+1:].strip()
ll=f.readline()
DateAndTime=ll[ll.index(":",)+1:].strip()
if args.debug:
    print("survey fund: %s"%SurveyName)
    print("Date of survey : %s"%DateAndTime)

ll=f.readline()
while len(ll)>0 and not ll.startswith("Electrodes Positions"):
        ll=f.readline()

if len(ll) <1:
    raise IOError("unable to find electrodes")
numStations=int(f.readline())

stations={}
for i in range(numStations):
    ll=f.readline()
    d=ll.split()
    sid=int(d[0])
    sx=float(d[1])
    sy=float(d[2])
    sz=float(d[3])
    stations[sid]=np.array([sx,sy,sz])

DistanceElectrodes=1e999
xmin=1e999
xmax=-1e999
ymin=1e999
ymax=-1e999
for i in stations:
    for j in stations:
        if i < j:
            d=np.linalg.norm(stations[j]-stations[i])
            DistanceElectrodes=min(DistanceElectrodes, d)
        xmin=min(xmin,stations[i][0])
        xmax=max(xmax,stations[i][0])
        ymin=min(ymin,stations[i][1])
        ymax=max(ymax,stations[i][1])

Core_x = xmax-xmin+2*max(xmax-xmin,ymax-ymin)*args.relRange/100.
Core_y = ymax-ymin+2*max(xmax-xmin,ymax-ymin)*args.relRange/100.
Core_z = max(xmax-xmin,ymax-ymin)*args.relRange/100.
if args.centralize:
    offsetX=-(xmax+xmin)/2.
    offsetY=-(ymax+ymin)/2.
    offsetZ=0
else:
    offsetX=0
    offsetY=0
    offsetZ=0

if args.filename is None:
    STFN=SurveyName+".loc"
else:
    STFN=args.stationfile 


if len(stations)>0:
    fs=open(STFN,'w')
    for sid in stations:
        fs.write("%s, %e, %e, %e\n"%(sid, offsetX+stations[sid][0], offsetY+stations[sid][1], offsetZ+stations[sid][2]))
if args.debug:
    print("%s stations found."%numStations)
    print("minimal distance of electrodes="%DistanceElectrodes)
    print("x range %s - %s"%(xmin+offsetX,xmax+offsetX))
    print("y range %s - %s"%(ymin+offsetY,ymax+offsetY))
    if len(stations)>0:
        print("stations are written to file %s."%STFN)



ll=f.readline()
while len(ll)>0 and not ll.startswith("Number of current pairs"):
        ll=f.readline()
if len(ll) <1:
    raise IOError("unable to find survey data")
numExperiments=int(f.readline())

ll=f.readline()
ll=f.readline()
ll=f.readline()

if args.debug:
    print("%s experiments found"%numExperiments)

## column 1 = dipole number
## column 3 = positive node
## column 2 = negative node
## column 4 = current
## columns 5 to  = 68 voltages at nodes 1, 3, 5,...63, 2, 4, ... 64 

WRAP=int(numStations/2)
index=[ WRAP* (g%2) + int(g/2)  for g in range(numStations) ]
#print index
records={}    
for i in range(numExperiments):
    ll=f.readline().split()
    ide=int(ll[0])
    A=int(ll[1])
    B=int(ll[2])
    I=float(ll[3])
    obs=ll[4:]
    if abs(I)>0:
        for N0 in range(numStations):
            for M0 in range(N0):
                #print ide, A, B, M0+1, N0+1, I, obs[index[A-1]], obs[index[B-1]], obs[index[M0]], obs[index[N0]]
                #if (A,B) ==(1, 64):
                #    print obs[index[M0]], obs[index[M0]]
                if not ( obs[index[M0]] == args.NaN or obs[index[N0]] == args.NaN ):
                    D=(float(obs[index[N0]])-float(obs[index[M0]]))/I
                    if abs(D)>0:
                        records[(A,B, M0+1, N0+1)]=D
                        
if args.filename is None:
    FN=SurveyName+".csv"
else:
    FN=args.filename 

if args.wenner:
    survey={}
    for a in range(1,numStations): # distance of central measuring electrodes
        for n in range(1,(numStations-a)/(2*a)+1): # n*a = distance of central measuring from charging electrodes
            for A in range(0,numStations-(2*n+1)*a): # first charging electrode
                M=A+n*a
                N=M+a
                B=N+n*a
                #print a,n, A," :", (A,M,N,B)
                if not B<numStations:
                    raise ValueError
                #if records.has_key((A+1,B+1,M+1,N+1)):
                if records.has_key((A+1,B+1,M+1,N+1)):
                    if records[(A+1,B+1,M+1,N+1)] <0 :
                        survey[(A+1,B+1,N+1, M+1)]=-records[(A+1,B+1,M+1,N+1)]
                    else:
                        survey[(A+1,B+1,M+1,N+1)]= records[(A+1,B+1,M+1,N+1)]
                elif records.has_key((A+1,B+1, N+1, M+1)):
                    if  records[(A+1,B+1,N+1, M+1)] <0 :
                        survey[(A+1,B+1,M+1,N+1)]= - records[(A+1,B+1,N+1, M+1)]
                    else:
                        survey[(A+1,B+1,N+1, M+1)] = records[(A+1,B+1,N+1, M+1)]
                elif records.has_key((B+1,A+1,M+1,N+1)):
                    if records[(B+1,A+1,M+1,N+1)] < 0:
                        survey[(A+1,B+1,M+1,N+1)]= - records[(B+1,A+1,M+1,N+1)]
                    else:
                        survey[(A+1,B+1,N+1, M+1)]= records[(B+1,A+1,M+1,N+1)]
                elif records.has_key((B+1,A+1, N+1, M+1)):
                    if records[(B+1,A+1,N+1, M+1)] < 0:
                        survey[(A+1,B+1,M+1,N+1)] =  records[(B+1,A+1,N+1, M+1)]
                    else:
                        survey[(A+1,B+1,N+1,M+1)] = -records[(B+1,A+1,N+1, M+1)]
else:
    survey=records
    
fs=open(FN,'w')
for A, B, M, N in survey:
    fs.write("%d, %d, %d, %d, %e\n"%(A, B, M, N, survey[(A,B,M,N)]))
if args.debug:
    print("%s measurements written to %s."%(len(survey), FN))
datafile=FN    

    

if not args.nogeo:
    if args.geofile is None:
        FN=SurveyName+".geo"
    else:
        FN=args.geofile 
    f=open(FN,'w')
    f.write("// Dimensions\n")
    f.write("surveyXmin = %g;\n"%xmin)
    f.write("surveyXmax = %g;\n"%xmax)
    f.write("surveyYmin = %g;\n"%ymin)
    f.write("surveyYmax = %g;\n"%ymax)
    f.write("// extra space to create core region:\n")
    f.write("surveyRangeX = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relRange/100.))
    f.write("surveyRangeY = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relRange/100.))
    f.write("surveyRangeZ = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relRange/100.))    
    f.write("// extra space to create core region:\n")
    f.write("surveyOffsetX = %g;\n"%(offsetX))
    f.write("surveyOffsetY = %g;\n"%(offsetY))
    f.write("surveyOffsetZ = %g;\n"%(offsetZ))    
    f.write("// extra space to for bounding box:\n")
    f.write("bbX = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relPadding/100.))
    f.write("bbY = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relPadding/100.))
    f.write("bbZ = %g;\n"%(max(xmax-xmin,ymax-ymin)*args.relPadding/100.))


    f.write("// element sizes\n")
    f.write("meshsizeBB = %g;\n"%(0.1*max(xmax-xmin+2*max(xmax-xmin,ymax-ymin)*args.relRange/100., ymax-ymin+2*max(xmax-xmin,ymax-ymin)*args.relRange/100.)))
    f.write("meshsizeCore = %g;\n"%DistanceElectrodes)
    f.write("meshsizeAtElectrodes = %g/2;\n"%DistanceElectrodes)
    

    
    f.write("\n")
    f.write("// core region:\n")
    f.write("Bx = surveyXmax-surveyXmin+2*surveyRangeX;\n")
    f.write("By = surveyYmax-surveyYmin+2*surveyRangeY;\n")
    f.write("Bz = surveyRangeZ;\n")
    f.write("bx = surveyXmin-surveyRangeX+surveyOffsetX;\n")
    f.write("by = surveyYmin-surveyRangeY+surveyOffsetY;\n")
    f.write("bz = -surveyRangeZ+surveyOffsetZ;\n")
    f.write("\n")
    f.write("// bounding box\n") 
    f.write("Wx = Bx+2*bbX;\n")
    f.write("Wy = By+2*bbY;\n")
    f.write("Wz = Bz+bbZ;\n")
    f.write("wx = bx-bbX;\n")
    f.write("wy = by-bbY;\n")
    f.write("wz = bz-bbZ;\n")

    f.write("\n")
    f.write("// Big Box top points\n")
    f.write("Point(1) = {  wx,    wy   , wz+Wz, meshsizeBB};\n")
    f.write("Point(2) = {  wx+Wx, wy   , wz+Wz, meshsizeBB};\n")
    f.write("Point(3) = {  wx,    wy+Wy,  wz+Wz, meshsizeBB};\n")
    f.write("Point(4) = {  wx+Wx, wy+Wy, wz+Wz, meshsizeBB};\n")
    f.write("\n")
    f.write("// Big Box bottom points\n")
    f.write("Point(5) = {  wx,    wy, wz, meshsizeBB};\n")
    f.write("Point(6) = {  wx+Wx, wy, wz, meshsizeBB};\n")
    f.write("Point(7) = {  wx,    wy+Wy, wz, meshsizeBB};\n")
    f.write("Point(8) = {  wx+Wx, wy+Wy, wz, meshsizeBB};\n")
    f.write("\n")
    f.write("// little box top points\n")
    f.write("Point(9)  = {  bx,     by   , bz+Bz, meshsizeCore};\n")
    f.write("Point(10) = {  bx+Bx,  by   , bz+Bz, meshsizeCore};\n")
    f.write("Point(11) = {  bx,     by+By, bz+Bz, meshsizeCore};\n")
    f.write("Point(12) = {  bx+Bx,  by+By, bz+Bz, meshsizeCore};\n")
    f.write("\n")
    f.write("// little box bottom points\n")
    f.write("Point(13) = {  bx,     by,    bz, meshsizeCore};\n")
    f.write("Point(14) = {  bx+Bx,  by,    bz, meshsizeCore};\n")
    f.write("Point(15) = {  bx,     by+By, bz, meshsizeCore};\n")
    f.write("Point(16) = {  bx+Bx,  by+By, bz, meshsizeCore};\n")
    f.write("\n")
    f.write("// all lines in positive axes directions\n")
    f.write("// Big Box top\n")
    f.write("Line(1) = {1,2} ;\n")
    f.write("Line(2) = {3,4} ;\n")
    f.write("Line(3) = {1,3} ;\n")
    f.write("Line(4) = {2,4} ;\n")
    f.write("// Big Box bottom\n")
    f.write("Line(5) = {5,6} ;\n")
    f.write("Line(6) = {7,8} ;\n")
    f.write("Line(7) = {5,7} ;\n")
    f.write("Line(8) = {6,8} ;\n")
    f.write("// Big Box sides\n")
    f.write("Line(9) = {1,5} ;\n")
    f.write("Line(10) = {2,6} ;\n")
    f.write("Line(11) = {3,7} ;\n")
    f.write("Line(12) = {4,8} ;\n")
    f.write("\n")
    f.write("// little box top \n")
    f.write("Line(13) = {9,10} ;\n")
    f.write("Line(14) = {11,12} ;\n")
    f.write("Line(15) = {9,11} ;\n")
    f.write("Line(16) = {10,12} ;\n")
    f.write("// little box bottom\n")
    f.write("Line(17) = {13,14} ;\n")
    f.write("Line(18) = {15,16} ;\n")
    f.write("Line(19) = {13,15} ;\n")
    f.write("Line(20) = {14,16} ;\n")
    f.write("// little box sides\n")
    f.write("Line(21) = {9,13} ;\n")
    f.write("Line(22) = {10,14} ;\n")
    f.write("Line(23) = {11,15} ;\n")
    f.write("Line(24) = {12,16} ;\n")
    f.write("\n")
    f.write("// Big Box surfaces minus little box at the top\n")
    f.write("// outward normals!!!\n")
    f.write("Line Loop(25) = {1,4,-2,-3} ;   // top    z = 0\n")
    f.write("Line Loop(31) = {13,16,-14,-15} ;   // top    z = 0\n")
    f.write("Line Loop(26) = {7,6,-8,-5} ;   // bottom z = WZ\n")
    f.write("Line Loop(27) = {9,5,-10,-1} ;  // front  y = -wy\n")
    f.write("Line Loop(28) = {2,12,-6,-11} ; // back   y =  wy\n")
    f.write("Line Loop(29) = {3,11,-7,-9} ;  // left   x = -wx\n")
    f.write("Line Loop(30) = {10,8,-12,-4} ; // right  x =  wx\n")
    f.write("\n")
    f.write("Plane Surface(1) = {25, -31};  // top    z = 0\n")
    f.write("Plane Surface(2) = {26};  // bottom z = WZ\n")
    f.write("Plane Surface(3) = {27};  // front  y = -wy\n")
    f.write("Plane Surface(4) = {28};  // back   y =  wy\n")
    f.write("Plane Surface(5) = {29};  // left   x = -wx\n")
    f.write("Plane Surface(6) = {30};  // right  x =  wx\n")
    f.write("\n")
    f.write("// little box surfaces \n")
    f.write("// outward normals!!!\n")
    f.write("Line Loop(32) = {19,18,-20,-17} ;   // bottom z = bz\n")
    f.write("Line Loop(33) = {21,17,-22,-13} ;   // front  y = -by\n")
    f.write("Line Loop(34) = {14,24,-18,-23} ;   // back   y =  by\n")
    f.write("Line Loop(35) = {15,23,-19,-21} ;   // left   x = -bx\n")
    f.write("Line Loop(36) = {22,20,-24,-16} ;   // right  x =  bx\n")
    f.write("\n")
    f.write("Plane Surface(7) = {31};   // top    z = 0\n")
    f.write("Plane Surface(8) = {32};   // bottom z = bz\n")
    f.write("Plane Surface(9) = {33};   // front  y = -by\n")
    f.write("Plane Surface(10) = {34};  // back   y =  by\n")
    f.write("Plane Surface(11) = {35};  // left   x = -bx\n")
    f.write("Plane Surface(12) = {36};  // right  x =  bx\n")
    f.write("\n")
    f.write("k = newp;\n")
    f.write("// electrodes  (in the little box)\n")

    n=0
    for e in stations:
            n+=1
            f.write("Point(k+%d) = { %e, %e, %e, meshsizeAtElectrodes};\n"%(n, offsetX+stations[e][0], offsetY+stations[e][1], offsetZ+stations[e][2]))
            f.write("Point{k+%d} In Surface{7};\n"%n)
    #        f.write("Physical Point(\"e%d\"") = { k+%d };\n"%(e,n))

    f.write("// Surface loops and volumes\n")
    f.write("\n")
    f.write("Surface Loop(1) = {1,2,3,4,5,6,8,9,10,11,12};\n")
    f.write("Surface Loop(2) = {7,8,9,10,11,12};\n")
    f.write("\n")
    f.write("Volume(1) = {1};\n")
    f.write("Volume(2) = {2};\n")
    f.write("\n")
    f.write("Physical Volume(\"padding\")  = { 1 } ;\n")
    f.write("Physical Surface(\"top\")    = { 1, 7 };\n")
    f.write("Physical Surface(\"bottom\") = { 2 };\n")
    f.write("Physical Surface(\"front\")  = { 3 };\n")
    f.write("Physical Surface(\"back\")   = { 4 };\n")

    f.write("Physical Surface(\"left\")  = { 5 };\n")
    f.write("Physical Surface(\"right\")   = { 6 };\n")

    f.write("Physical Volume(\"core_region\")  = { 2 } ;\n")
    f.close()
    #f.write("Physical Surface("lbTop")    = { 7 };
    #f.write("Physical Surface("lbBottom") = { 8 };
    #f.write("Physical Surface("lbFront")  = { 9 };
    #f.write("Physical Surface("lbBack")   = { 10 };
    #f.write("Physical Surface("lbLeft")  = { 11 };
    #f.write("Physical Surface("lbRight")   = { 12 };
    print "run 'gmsh -3 -optimize -o "+SurveyName+".msh "+FN+"'"



if not args.configfile:
    if args.configfile is None:
        FN=SurveyName+".Ypadding"
    else:
        FN=args.configfile 
    f=open(FN,'w')
    f.write("PROJECTNAME='%s'\n"%SurveyName)
    f.write("meshfile='%s'\n"%(SurveyName+".fly"))
    f.write("datafile='%s'\n"%datafile)
    f.write("datadelimiter=','\n")
    f.write("stationfile='%s'\n"%(STFN))
    f.write("stationdelimiter=','\n")
    f.write("stationsFMT='%s'\n")
    
    f.write("# Size of the core region\n")
    f.write("Core_x = %g\n"%Core_x)
    f.write("Core_y = %g\n"%Core_y)
    f.write("Core_z = %g\n"%Core_z)

    f.write("# grid for defining conductivity in the core region:\n")
    f.write("Ncell_x=400\n")
    f.write("Ncell_z=200\n")
    f.write("sigma_BB=0.001\n")
    f.write("sigma_Core=0.01\n")

    f.write("# electrodes on surface z=0\n")
    f.write("# need to be within [-Core_x/2, Core_x] X  [-Core_y/2, Core_y]\n")
    f.write("# dictionary with { <id>: <location>}\n")
    f.write("electrodes = {  \n")
    nn=0
    for e in stations:
        f.write( " %d : (%g, %g)"%(e, offsetX+stations[e][0], offsetY+stations[e][1]))
        if nn < len(stations)-1:
            f.write( ",\n")
        else:
            f.write( "}\n")
        nn+=1
    f.write("numElectrodes = len(electrodes)\n")

    f.write("usesStationCoordinates=False\n") 
    f.write("useLogDefect=True\n")

    f.write("dipoleInjections=True\n")
    f.write("dipoleMeasurements=True\n")
    f.write("region_fixed=['padding']\n")

    f.write("sigma_background=0.001\n")
    f.write("sigma0={'core_region': sigma_Core, 'padding' : sigma_BB}\n")
    f.write("L_stations=%g\n"%(DistanceElectrodes/5))

    f.write("outfile='sigma'\n")
    f.write("alpha0=1.\n")
    f.write("alpha1=(%g)**2/2\n"%DistanceElectrodes)
    f.write("w0=0.\n")
    f.write("w1=1e-3\n")
    f.close()
    if args.debug:
        print("config file %s generated."%(FN))
