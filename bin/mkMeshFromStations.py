#!/usr/bin/python3
"""
fingal - Creates a mesh fly file using the station location information. 

"""
import argparse
import importlib, sys, os
from fingal import *
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
sys.path.append(os.getcwd())
#sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='Creates a mesh fly file using the station location information. the gmsh mesh generator is used.', epilog="fingal by l.gross@uq.edu.auversion 21/12/2020")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--coredepth', '-d',  dest='coredepth', type=int, default=45, help="core depth relative to core width around core in %% of electrode area size (default 45)")
parser.add_argument('--extracore', '-e',  dest='extracore', type=int, default=20, help="relative extracore padding within the core around electrodes in %%  of edge legth (default 20)")
parser.add_argument('--padding', '-p',  dest='padding', type=int, default=150, help="relative padding around core in %% of core eddge length (default 150)")
parser.add_argument('--coremeshfactor', '-C',  dest='coremeshfactor', type=float, default=1., help="refinement factor of mesh in core relative to the electrode distance (default 0.1)")
parser.add_argument('--stationmeshfactor', '-s',  dest='stationmeshfactor', type=float, default=0.2, help="refinement factor at stations relative to core mesh size (default 0.2)")
parser.add_argument('--paddingmesh', '-P',  dest='paddingmesh', type=float, default=20, help="number of element on the longest edge of the padding region (default 20)")
parser.add_argument('--geo', '-g',  dest='geofile', type=str, default="tmp", help="name of gmsh geofile to generate")
parser.add_argument('--noMesh', '-N',  dest='mshno', action='store_true', default=False, help="if set only gmsh geo file is generated but mesh generation is not started. Useful for debugging.")
parser.add_argument('--silofile', '-o', dest='silofile', metavar='SILOFILE', type=str, default='mesh', help='name of silo output file for mesh (generated)')

args = parser.parse_args()

print("** This generates a mesh fly file for inversion from a station file **")

config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

positions=np.array([ [ elocations[s][0], elocations[s][1], elocations[s][2]] for s in elocations ] )
zCore=np.mean(positions[:,2])
zmin=positions[:,2].min()
xmin=positions[:,0].min()
xmax=positions[:,0].max()
ymin=positions[:,1].min()
ymax=positions[:,1].max()

DistanceElectrodes=1e999
for i in range(len(elocations)):
    for j in range(len(elocations)):
        if i < j:
            d=np.linalg.norm(positions[j]-positions[i])
            DistanceElectrodes=min(DistanceElectrodes, d)

print("level of survey = ", zCore)
print("x range = ", xmin, xmax)
print("y range = ", ymin, ymax)
print("DistanceElectrodes =",DistanceElectrodes)
diagonalAreaOfElectrodes= ((xmax - xmin) ** 2 + (ymax - ymin) ** 2) ** 0.5
assert diagonalAreaOfElectrodes > 0., "area of electrodes is zero."
XextensionCore = (max(xmax - xmin, ymax - ymin) * args.extracore) / 100.
YextensionCore = XextensionCore

if xmax-xmin < 1e-5 * diagonalAreaOfElectrodes:
  XextensionCore=YextensionCore
if ymax-ymin < 1e-5 * diagonalAreaOfElectrodes:
  YextensionCore=XextensionCore

CoreThickness= (diagonalAreaOfElectrodes * args.coredepth) / 100.
#zmincore=min(zmin, (fit[0]+fit[1]*(xmin-XextensionCore)+fit[2]*(ymin-YextensionCore)), fit[0]+fit[1]*(xmax+XextensionCore)+fit[2]*(ymin-YextensionCore), fit[0]+fit[1]*(xmax+XextensionCore)+fit[2]*(ymax+YextensionCore), fit[0]+fit[1]*(xmin-XextensionCore)+fit[2]*(ymax+YextensionCore))

Xpadding= ((xmax - xmin + 2 * XextensionCore) * args.padding) / 100.
Ypadding= ((ymax - ymin + 2 * YextensionCore) * args.padding) / 100.
Zpadding=min(Xpadding, Ypadding)  
Xpadding,Ypadding=Zpadding,Zpadding

XminBB = xmin - XextensionCore - Xpadding
XmaxBB = xmax + XextensionCore + Xpadding
YminBB = ymin - YextensionCore - Ypadding
YmaxBB = ymax + YextensionCore + Ypadding
ZminBB = zCore - CoreThickness - Zpadding
ZminCore = zCore - CoreThickness


out=""
out+="Mesh.MshFileVersion = 2.2;\n"
out+="// Core:\n"
out+="XminCore = %s;\n"%(xmin - XextensionCore)
out+="XmaxCore = %s;\n"%(xmax + XextensionCore)
out+="YminCore = %s;\n"%(ymin - YextensionCore)
out+="YmaxCore = %s;\n"%(ymax + YextensionCore)
out+="ZCore = %s;\n"%zCore
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
#out+="meshSizeCore = %s/%s;\n"%(max(xmax - xmin + 2 * XextensionCore, ymax - ymin + 2 * YextensionCore), args.coremesh)
#out+="meshSizeBB = %s/%s;\n"%(max(xmax - xmin + 2 * XextensionCore + 2 * Xpadding, ymax - ymin + 2 * YextensionCore + 2 * Ypadding), args.paddingmesh)
#out+="meshSizeElectrodes = meshSizeCore*%s;\n"%(args.stationmeshfactor)
######
out += "DistanceElectrodes = %g;\n" % (DistanceElectrodes)
out += "meshSizeCore = %g * DistanceElectrodes;\n" % (args.coremeshfactor)
out += "meshSizeBB = %g/%d;\n" % (
max(xmax - xmin + 2 * XextensionCore + 2 * Xpadding, ymax - ymin + 2 * YextensionCore + 2 * Ypadding),
args.paddingmesh)
out += "meshSizeElectrodes = DistanceElectrodes * %s;\n" % (args.stationmeshfactor)
###

out+="\n"
out+="Point(1) = {XminCore, YminCore, ZCore, meshSizeCore};\n"
out+="Point(2) = {XmaxCore, YminCore, ZCore, meshSizeCore};\n"
out+="Point(3) = {XmaxCore, YmaxCore, ZCore, meshSizeCore};\n"
out+="Point(4) = {XminCore, YmaxCore, ZCore, meshSizeCore};\n"
out+="\n"
out+="Point(5) = {XminCore, YminCore, ZminCore, meshSizeCore};\n"
out+="Point(6) = {XmaxCore, YminCore, ZminCore, meshSizeCore};\n"
out+="Point(7) = {XmaxCore, YmaxCore, ZminCore, meshSizeCore};\n"
out+="Point(8) = {XminCore, YmaxCore, ZminCore, meshSizeCore};\n"
out+="\n"
out+="Point(9) = {XminBB, YmaxBB, ZCore, meshSizeBB };\n"
out+="Point(10) = {XminBB, YminBB, ZCore, meshSizeBB };\n"
out+="Point(11) = {XmaxBB, YminBB, ZCore, meshSizeBB };\n"
out+="Point(12) = {XmaxBB, YmaxBB, ZCore, meshSizeBB };\n"
out+="\n"
out+="Point(17) = {XminBB, YmaxBB, ZminBB, meshSizeBB };\n"
out+="Point(18) = {XminBB, YminBB, ZminBB, meshSizeBB };\n"
out+="Point(19) = {XmaxBB, YminBB, ZminBB, meshSizeBB };\n"
out+="Point(20) = {XmaxBB, YmaxBB, ZminBB, meshSizeBB };\n"
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
out+="Plane Surface(7) = {-7};\n"
out+="Line Loop(8) = {19, 3, -21, -20};\n"
out+="Plane Surface(8) = {-8};\n"
out+="Line Loop(9) = {23, 21, 4, 12};\n"
out+="Plane Surface(9) = {-9};\n"
out+="Line Loop(10) = {12, 22, -11, -5};\n"
out+="Plane Surface(10) = {10};\n"
out+="Line Loop(11) = {1, -5, -4, -3};\n"
out+="Line Loop(12) = {9, 6, 7, 8};\n"
out+="Plane Surface(11) = {12};\n"
out+="Plane Surface(12) = {11, 12};\n"

out+="// electrodes  (in the core)\n"
out+="k=newp;\n"
for i,s in enumerate(elocations):
  out+="Point(k+%s)={ %s, %s, ZCore, meshSizeElectrodes};\n"%(i+1, elocations[s][0], elocations[s][1])
  out+="Point{k+%s} In Surface{11};\n"%(i+1)  
  out+='Physical Point("s%s")  = { k+%s } ;\n'%(s,i+1)
out+='Physical Surface("'+ config.faces[0] + '") = {6, 7, 8, 9, 10};\n'


out+="Surface Loop(1) = {1,2,3,4,5,11};\n"
out+="Volume(1) = {-1};\n"
out+="Surface Loop(2) = {6, -8, -9, -7, 10, 2, 1, 12, 5, 4, 3};\n"
out+="Volume(2) = {2};\n"
out+='Physical Volume("padding")  = { 2 } ;\n'
out+='Physical Volume("'+ config.core[0] + '")  = { 1 } ;\n'

GEOFN2=args.geofile+".geo"
MSHN3=args.geofile+".msh"


open(GEOFN2,'w').write(out)
print("3D geometry has been written to %s"%GEOFN2)

if not args.mshno:
    rp=subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHN3, GEOFN2])
    #rp.check_returncode()
    print(">> GMSH mesh file %s was generated."%MSHN3)


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
    if args.silofile:
        saveSilo(args.silofile , tag = makeTagMap(Function(domain)))
        print(f'Mesh written to file {args.silofile}.silo.')
    
