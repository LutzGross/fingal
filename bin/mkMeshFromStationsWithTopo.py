#!/usr/bin/python3
"""
fingal - Creates a mesh fly file using the station location information. 

"""
import argparse
import importlib, sys, os
from fingal import *
import subprocess
from esys.finley import ReadGmsh
sys.path.append(os.getcwd())
#sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='Creates a mesh fly file using the station location information. the gmsh mesh generator is used.', epilog="fingal by l.gross@uq.edu.auversion 21/12/2020")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--topo', '-t',  dest='topo', action='store_true', default=False, help="topography is added from station locations (experimental)")
parser.add_argument('--inner', '-i',  dest='inner', type=int, default=20, help="relative inner padding within the core around electrodes in %% (default 20)")
parser.add_argument('--outer', '-o',  dest='outer', type=int, default=80, help="relative outer padding around core in %% (default 80)")
parser.add_argument('--depth', '-d',  dest='depth', type=int, default=45, help="depth relative to core width around core in %% (default 45)")
parser.add_argument('--plot', '-P',  dest='plot', type=str, default=None, help="file name to plot topography")
parser.add_argument('--coremesh', '-c',  dest='coremesh', type=float, default=50, help="number of element on the longest edge of the core region (default 50)")
parser.add_argument('--stationmesh', '-s',  dest='stationmesh', type=float, default=0.2, help="refinement factor at stations relative to core mesh size (default 0.3)")
parser.add_argument('--paddingmesh', '-p',  dest='paddingmesh', type=float, default=20, help="number of element on the longest edge of the padding region (default 20)")
parser.add_argument('--geo', '-g',  dest='geofile', type=str, default="tmp", help="name of gmsh geofile to generate")
parser.add_argument('--flyno', '-f',  dest='flyno', action='store_true', default=False, help="if set no fly file conversion is started.")
parser.add_argument('--mshno', '-m',  dest='mshno', action='store_true', default=False, help="if set only gmsh geo file is generated but mesh generation is not started. Useful for debugging.")

args = parser.parse_args()

print("** This generates a mesh fly file **")
if args.topo:
    print("Topography is added.")


config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

positions=np.array([ [ elocations[s][0], elocations[s][1], elocations[s][2]] for s in elocations ] )
zpos=np.mean(positions[:,2])
zmin=positions[:,2].min()
xmin=positions[:,0].min()
xmax=positions[:,0].max()
ymin=positions[:,1].min()
ymax=positions[:,1].max()
fit , res, r, s =np.linalg.lstsq(np.stack((np.ones((len(positions[:,0]),)), positions[:,0], positions[:,1])).T, positions[:,2], rcond=None)
#fit[0]=zpos
fit[1]=0
fit[2]=0
fit[0]=zpos
#print(fit)
print("mean of topograpy = ", zpos)
print("x range = ", xmin, xmax)
print("y range = ", ymin, ymax)

ebox=max(xmax-xmin, ymax-ymin)
assert ebox > 0., "area of electrodes is zero."
fx=((xmax-xmin)*args.inner)/100.
fy=((ymax-ymin)*args.inner)/100.

if xmax-xmin < 1e-5 * ebox:
  fx=fy  
if ymax-ymin < 1e-5 * ebox:
  fy=fx  


fz=(min(xmax-xmin, ymax-ymin)*args.depth)/100.
zmincore=min(zmin, (fit[0]+fit[1]*(xmin-fx)+fit[2]*(ymin-fy)), fit[0]+fit[1]*(xmax+fx)+fit[2]*(ymin-fy), fit[0]+fit[1]*(xmax+fx)+fit[2]*(ymax+fy), fit[0]+fit[1]*(xmin-fx)+fit[2]*(ymax+fy))
px=((xmax-xmin+2*fx)*args.outer)/100.
py=((ymax-ymin+2*fy)*args.outer)/100.
pz=min(px, py)  
px,py=pz,pz

XminBB = xmin-fx-px
XmaxBB = xmax+fx+px
YminBB = ymin-fy-py
YmaxBB = ymax+fy+py
ZminBB = zmincore-fz-pz
ZminCore = (zmincore-fz)

out=""
out+="Mesh.MshFileVersion = 2.2\n"
out+="// Core:\n"
out+="XminCore = %s;\n"%(xmin-fx)
out+="XmaxCore = %s;\n"%(xmax+fx)
out+="YminCore = %s;\n"%(ymin-fy)
out+="YmaxCore = %s;\n"%(ymax+fy)
out+="ZminCore = %s;\n"%ZminCore
out+="\n"
out+="// big box\n"
out+="XminBB = %s;\n"%XminBB
out+="XmaxBB = %s;\n"%XmaxBB
out+="YminBB = %s;\n"%YminBB
out+="YmaxBB = %s;\n"%YmaxBB
out+="ZminBB = %s;\n"%ZminBB
out+="\n"

out+="// element sizes\n"
out+="mshC = %s/%s;\n"%(max(xmax-xmin+2*fx,ymax-ymin+2*fy), args.coremesh)
out+="mshBB = %s/%s;\n"%(max(xmax-xmin+2*fx+2*px,ymax-ymin+2*fy+2*py), args.paddingmesh)
out+="mshE = mshC*%s;\n"%(args.stationmesh)
out+="\n"
out+="Point(1) = {XminCore, YminCore, %s, mshC};\n"%(fit[0]+fit[1]*(xmin-fx)+fit[2]*(ymin-fy))
out+="Point(2) = {XmaxCore, YminCore, %s, mshC};\n"%(fit[0]+fit[1]*(xmax+fx)+fit[2]*(ymin-fy))
out+="Point(3) = {XmaxCore, YmaxCore, %s, mshC};\n"%(fit[0]+fit[1]*(xmax+fx)+fit[2]*(ymax+fy))
out+="Point(4) = {XminCore, YmaxCore, %s, mshC};\n"%(fit[0]+fit[1]*(xmin-fx)+fit[2]*(ymax+fy))
out+="\n"
out+="Point(5) = {XminCore, YminCore, ZminCore, mshC};\n"
out+="Point(6) = {XmaxCore, YminCore, ZminCore, mshC};\n"
out+="Point(7) = {XmaxCore, YmaxCore, ZminCore, mshC};\n"
out+="Point(8) = {XminCore, YmaxCore, ZminCore, mshC};\n"
out+="\n"
out+="Point(9) = {XminBB, YmaxBB, %s, mshBB };\n"%(fit[0]+fit[1]*XminBB+fit[2]*YmaxBB)
out+="Point(10) = {XminBB, YminBB, %s, mshBB };\n"%(fit[0]+fit[1]*XminBB+fit[2]*YminBB)
out+="Point(11) = {XmaxBB, YminBB, %s, mshBB };\n"%(fit[0]+fit[1]*XmaxBB+fit[2]*YminBB)
out+="Point(12) = {XmaxBB, YmaxBB, %s, mshBB };\n"%(fit[0]+fit[1]*XmaxBB+fit[2]*YmaxBB)
out+="\n"
out+="Point(17) = {XminBB, YmaxBB, ZminBB, mshBB };\n"
out+="Point(18) = {XminBB, YminBB, ZminBB, mshBB };\n"
out+="Point(19) = {XmaxBB, YminBB, ZminBB, mshBB };\n"
out+="Point(20) = {XmaxBB, YmaxBB, ZminBB, mshBB };\n"
out+="\n"
out+="Line(1) = {10, 9};\n"
out+="Line(3) = {10, 11};\n"
out+="Line(4) = {11, 12};\n"
out+="Line(5) = {12, 9};\n"
out+="Line(6) = {4, 3};\n"
out+="Line(7) = {3, 2};\n"
out+="Line(8) = {2, 1};\n"
out+="Line(9) = {1, 4};\n"
out+="Line(10) = {4, 8};\n"
out+="Line(11) = {9, 17};\n"
out+="Line(12) = {12, 20};\n"
out+="Line(13) = {3, 7};\n"
out+="Line(14) = {2, 6};\n"
out+="Line(15) = {6, 5};\n"
out+="Line(16) = {1, 5};\n"
out+="Line(17) = {5, 8};\n"
out+="Line(18) = {17, 18};\n"
out+="Line(19) = {18, 10};\n"
out+="Line(20) = {18, 19};\n"
out+="Line(21) = {19, 11};\n"
out+="Line(22) = {20, 17};\n"
out+="Line(23) = {20, 19};\n"
out+="Line(24) = {8, 7};\n"
out+="Line(25) = {7, 6};\n"
out+="Line Loop(1) = {6, 13, -24, -10};\n"
out+="Plane Surface(1) = {1};\n"
out+="Line Loop(2) = {10, -17, -16, 9};\n"
out+="Plane Surface(2) = {2};\n"
out+="Line Loop(3) = {15, 17, 24, 25};\n"
out+="Plane Surface(3) = {3};\n"
out+="Line Loop(4) = {13, 25, -14, -7};\n"
out+="Plane Surface(4) = {4};\n"
out+="Line Loop(5) = {14, 15, -16, -8};\n"
out+="Plane Surface(5) = {5};\n"
out+="Line Loop(6) = {19, 1, 11, 18};\n"
out+="Plane Surface(6) = {6};\n"
out+="Line Loop(7) = {22, 18, 20, -23};\n"
out+="Plane Surface(7) = {7};\n"
out+="Line Loop(8) = {19, 3, -21, -20};\n"
out+="Plane Surface(8) = {8};\n"
out+="Line Loop(9) = {23, 21, 4, 12};\n"
out+="Plane Surface(9) = {9};\n"
out+="Line Loop(10) = {12, 22, -11, -5};\n"
out+="Plane Surface(10) = {10};\n"
out+="Line Loop(11) = {1, -5, -4, -3};\n"
out+="Line Loop(12) = {9, 6, 7, 8};\n"
out+="Plane Surface(11) = {12};\n"
out+="Plane Surface(12) = {11, 12};\n"

out+="// electrodes  (in the core)\n"
out+="k=newp;\n"
for i,s in enumerate(elocations):
  out+="Point(k+%s)={ %s, %s, %s, mshE};\n"%(i+1, elocations[s][0], elocations[s][1], fit[0]+fit[1]*elocations[s][0]+fit[2]*elocations[s][1])
  out+="Point{k+%s} In Surface{11};\n"%(i+1)  
  out+='Physical Point("s%s")  = { k+%s } ;\n'%(s,i+1)
out+='Physical Surface("faces") = { 1, 2, 3, 4,5 ,6,7,8,9,10, 11,12};\n'

GEOFN2=args.geofile+".geo"
MSHN3=args.geofile+".msh"


if not args.topo:
    out+="Surface Loop(1) = {1,2,3,4,5,11};\n"
    out+="Volume(1) = {-1};\n"
    out+="Surface Loop(2) = {6, 8, 9, 7, 10, 2, 1, 12, 5, 4, 3};\n"
    out+="Volume(2) = {2};\n"
    out+='Physical Volume("padding")  = { 2 } ;\n'
    out+='Physical Volume("core")  = { 1 } ;\n'
    open(GEOFN2,'w').write(out)
    print("3D geometry has been written to %s"%GEOFN2)
    if not args.mshno:
        rp=subprocess.run(["gmsh", "-3",  "-algo", "auto", "-o", MSHN3, GEOFN2])
        rp.check_returncode()
        print(">> GMSH mesh file %s was generated."%MSHN3)
else:
    # this is still experimental!
    GEOFN1=args.geofile+"_2D.geo"
    MSHN1=args.geofile+"_2D.msh"
    MSHN2=args.geofile+"_2D_topo.msh"


    open(GEOFN1,'w').write(out)
    print("surface geometry was written to %s"%GEOFN1)


    out="// Core:\n"
    out+="XminCore = %s;\n"%(xmin-fx)
    out+="XmaxCore = %s;\n"%(xmax+fx)
    out+="YminCore = %s;\n"%(ymin-fy)
    out+="YmaxCore = %s;\n"%(ymax+fy)
    out+="ZminCore = %s;\n"%ZminCore
    out+="// element sizes\n"
    out+="mshC = %s/%s;\n"%(min(xmax-xmin+2*fx,ymax-ymin+2*fy), args.coremesh)
    out+="mshBB = %s/%s;\n"%(min(xmax-xmin+2*fx+2*px,ymax-ymin+2*fy+2*py), args.paddingmesh)

    out+='Merge "%s"; // merge modified msh\n'%MSHN2

    out+="\n"
    out+="Surface Loop(1) = {1,2,3,4,5,11};\n"
    out+="Volume(1) = {-1};\n"
    out+="Surface Loop(2) = {6, 8, 9, 7, 10, 2, 1, 12, 5, 4, 3};\n"
    out+="Volume(2) = {2};\n"
    out+='Physical Volume("padding")  = { 2 } ;\n'
    out+='Physical Volume("core")  = { 1 } ;\n'
    out+="\n"
    out+="Field[1] = Box;\n"
    out+="Field[1].VIn = mshC;\n"
    out+="Field[1].VOut = mshBB;\n"
    out+="Field[1].XMin = XminCore;\n"
    out+="Field[1].XMax = XmaxCore;\n"
    out+="Field[1].YMin = YminCore;\n"
    out+="Field[1].YMax = YmaxCore;\n"
    out+="Field[1].ZMin = ZminCore;\n"
    out+="Field[1].ZMax = 1e55;\n"
    out+="Background Field = 1;\n"
    out+="\n"
    out+="Mesh.CharacteristicLengthExtendFromBoundary = 0;\n"
    out+="Mesh.CharacteristicLengthFromPoints = 0;\n"
    out+="Mesh.CharacteristicLengthFromCurvature = 0;\n"
    out+="Mesh.CharacteristicLengthMin=0.;\n"
    out+="Mesh.CharacteristicLengthMax=1e33;\n"

    open(GEOFN2,'w').write(out)
    print("3D geometry was written to %s"%GEOFN2)

    print("generate surface mesh:")


    rp=subprocess.run(["gmsh", "-2", "-algo", "auto", "-o", MSHN1, GEOFN1])
    rp.check_returncode()
    print(">> 2D GMSH mesh file %s generated."%MSHN1)

    #   new the surface at z=zpos    
    from scipy.interpolate import LinearNDInterpolator as Interpolator
    #from scipy.interpolate import CloughTocher2DInterpolator as Interpolator



    # create an interpolationn table:
    points = np.vstack((positions[:,:2], [ [XminBB, YmaxBB], [XminBB, YminBB ] , [XmaxBB, YminBB], [XmaxBB, YmaxBB] ]  ) )
    elevation = np.vstack( (positions[:,2:], [ [zmin], [zmin], [zmin], [zmin] ]))
    interpolator = Interpolator(points, elevation)
    # a quick plot:
    if args.plot:
        X = np.linspace(XminBB, XmaxBB, num=200)
        Y = np.linspace(YminBB, YmaxBB, num=200)
        X, Y = np.meshgrid(X, Y)
        Z = interpolator(X, Y)
        import matplotlib.pyplot as plt
        plt.figure()


        plt.pcolormesh(X, Y, Z[:,:,0])

        plt.colorbar() # Color Bar
        plt.scatter(points[:,0], points[:,1], s=4)
        plt.savefig(args.plot)
        print("topography pic written to %s"%args.plot)

    # apply topography:
    fin=open(MSHN1, 'r')
    fout=open(MSHN2,'w')
    line=fin.readline()
    while line:
        if line.startswith("$NOD") or line.startswith("$Nodes"):
                process_nodes=True
                numNodes=int(fin.readline())
                if line.startswith("$NOD"):
                    fout.write("$NOD\n")
                else:
                    fout.write("$Nodes\n")
                fout.write("%d\n"%numNodes)
                print(numNodes," nodes found.")
                cc=0
                while cc < numNodes:
                    line=fin.readline()
                    i,x,y,z=line.split(" ")
                    level=fit[0]
                    if float(z) > ZminCore :
                        znew=(interpolator(float(x), float(y))[0]-ZminCore)*(float(z)-ZminCore)/(level-ZminCore)+ZminCore
                    else:
                        znew=float(z)
                    fout.write("%s %g %g %g\n"%(i, float(x), float(y), znew ))
                    cc+=1
        else:
                fout.write(line)
        line=fin.readline()

    fin.close()
    fout.close()
    print("topography added to file %s and written to file %s."%(MSHN1, MSHN2))
    # now we are ready to generate the 3D mesh:  
    rp=subprocess.run(["gmsh", "-3",  "-algo", "frontal", "-o", MSHN3, GEOFN2])
    rp.check_returncode()
    print(">> GMSH mesh file %s generated."%MSHN3)


if not args.flyno and not args.mshno:
    dts=[]
    dps=[]
    for s in elocations:
        if config.stationsFMT:
            dts.append(config.stationsFMT%s)
        else:
            dts.append(s)
        dps.append(elocations[s])
    domain=ReadGmsh(MSHN3, 3, diracPoints=dps, diracTags=dts, optimize=True )
    domain.write(config.meshfile)
    print("Mesh written to file %s"%(config.meshfile))


    
