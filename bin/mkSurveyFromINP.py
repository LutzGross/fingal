#!/usr/bin/python3
"""
fingal - Creates a mesh fly file using the station location information.

"""
import argparse
import importlib, sys, os
from fingal import *
from datetime import datetime
import subprocess
from esys.finley import ReadGmsh
from esys.escript import Function, makeTagMap
from esys.weipa import saveSilo
sys.path.append(os.getcwd())
# sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(
    description='Creates a mesh fly file, station location information and survey (gmsh is used).',
    epilog="fingal by l.gross@uq.edu.auversion 21/12/2020")
parser.add_argument(dest='INP', metavar='INPFile', type=str, help='INP file')
parser.add_argument('--topo', '-t', dest='topo', action='store_true', default=False,
                    help="topography is added from station locations (ignored)")
parser.add_argument('--core', '-i', dest='core', type=int, default=20,
                    help="core around area of electrodes in %% of its diameter (default 20)")
parser.add_argument('--paddingh', '-H', dest='paddingh', type=int, default=100,
                    help="relative horizontal padding padding around core in %% of core (default 100)")
parser.add_argument('--paddingz', '-Z', dest='paddingz', type=int, default=100,
                    help="relative vertical padding padding around core in %% of core (default 100)")
parser.add_argument('--depth', '-d', dest='depth', type=int, default=45,
                    help="depth relative to core width around core in %% (default 45)")

#parser.add_argument('--plot', '-P', dest='plot', type=str, default=None, help="file name to plot topography")
parser.add_argument('--coremeshfactor', '-c', dest='coremeshfactor', type=float, default=1.,
                    help="factor for element size relative electrode distance (default 1)")
parser.add_argument('--stationmeshfactor', '-s', dest='stationmeshfactor', type=float, default=0.4,
                    help="refinement factor at stations relative to electrode distance (default 0.3)")
parser.add_argument('--paddingmesh', '-p', dest='paddingmesh', type=float, default=10,
                    help="number of element on the longest edge of the padding region (default 10)")
parser.add_argument('--geo', '-G', dest='geofile', type=str, default="tmp", help="name of gmsh geofile to generate")
parser.add_argument('--mesh', '-M', dest='flyfile', type=str, default="mesh.fly", help="name of fly mesh file to be generated")
parser.add_argument('--data', '-D', dest='datafile', type=str, default="data.csv", help="name of data file to be generated")
parser.add_argument('--stations', '-S', dest='stationfile', type=str, default="stations.csv", help="name of station/electrode location file to be generated")
parser.add_argument('--config', '-C', dest='configfile', type=str, default="config", help="name of the configuration file to be generated.")

args = parser.parse_args()
CDTEMPLATE=os.path.dirname(__file__)
stationsFMT="s%s"
silofile="mesh"

print("** This generate ERT inversion input from INP file **")
if args.topo:
    print("Topography is added (ignored).")


# read electrodes:
INPUT=open(args.INP, 'r')
print("input file ", args.INP)
for i in range(7):
    INPUT.readline()
# read stations:
NumElectrodes=int(INPUT.readline())
print("Number of Electrodes found is ", NumElectrodes)
positions=[]
electrode_in_use={}
electrodes_id=[]

for i in range(NumElectrodes):
    s=INPUT.readline().split()
    id=int(s[0])
    electrodes_id.append(id)
    positions.append([float(s[1]), float(s[2]), float(s[3])])
    electrode_in_use[id] = s[4] == '0'

positions = np.array(positions)
zpos = np.mean(positions[:, 2])
zmin = positions[:, 2].min()
zmax = positions[:, 2].max()
xmin = positions[:, 0].min()
xmax = positions[:, 0].max()
ymin = positions[:, 1].min()
ymax = positions[:, 1].max()
print(f"Area of electrodes [{xmin},{xmax}] x [{ymin}, {ymax}]")

diameter_area_electrodes =np.sqrt((xmax - xmin)**2 +(ymax - ymin)**2)
assert diameter_area_electrodes > 0., "area of electrodes is zero."

distmin=diameter_area_electrodes
for i in range(NumElectrodes-1):
   distmin=min(distmin, np.linalg.norm(positions[i]-positions[i+1:],axis=1).min())
print(f"minimum electrode distance = {distmin}")

fx = ((xmax - xmin) * args.core) / 100.
fy = ((ymax - ymin) * args.core) / 100.
fz = (diameter_area_electrodes * args.depth) / 100.

if xmax - xmin < 1e-5 * diameter_area_electrodes:
   fx = fy
if ymax - ymin < 1e-5 * diameter_area_electrodes:
   fy = fx

px = ((xmax - xmin + 2 * fx) * args.paddingh) / 100.
py = ((ymax - ymin + 2 * fy) * args.paddingh) / 100.
pz = ((zmax - zmin + fz) * args.paddingz) / 100.

pp = min(px, py, pz)
px, py = pp, pp

XminBB = xmin - fx - px
XmaxBB = xmax + fx + px
YminBB = ymin - fy - py
YmaxBB = ymax + fy + py
ZminBB = zpos - fz - pz
ZminCore = zpos - fz

out = ""
out += "Mesh.MshFileVersion = 2.2;\n"
out += "// Core:\n"
out += "Zpos = %s;\n" % zpos
out += "XminCore = %s;\n" % (xmin - fx)
out += "XmaxCore = %s;\n" % (xmax + fx)
out += "YminCore = %s;\n" % (ymin - fy)
out += "YmaxCore = %s;\n" % (ymax + fy)
out += "ZminCore = %s;\n" % ZminCore
out += "\n"
out += "// big box\n"
out += "XminBB = %s;\n" % XminBB
out += "XmaxBB = %s;\n" % XmaxBB
out += "YminBB = %s;\n" % YminBB
out += "YmaxBB = %s;\n" % YmaxBB
out += "ZminBB = %s;\n" % ZminBB
out += "\n"

out += "// element sizes\n"
out += "mshC = %s * %s;\n" % (distmin, args.coremeshfactor)
out += "mshBB = %s/%s;\n" % (max(xmax - xmin + 2 * fx + 2 * px, ymax - ymin + 2 * fy + 2 * py, zmax - zmin + fz + pz), args.paddingmesh)
out += "mshE = %s * %s;\n" % (distmin, args.stationmeshfactor)
out += "\n"
out += "Point(1) = {XminCore, YminCore, Zpos, mshC};\n"
out += "Point(2) = {XmaxCore, YminCore, Zpos, mshC};\n"
out += "Point(3) = {XmaxCore, YmaxCore, Zpos, mshC};\n"
out += "Point(4) = {XminCore, YmaxCore, Zpos, mshC};\n"
out += "\n"
out += "Point(5) = {XminCore, YminCore, ZminCore, mshC};\n"
out += "Point(6) = {XmaxCore, YminCore, ZminCore, mshC};\n"
out += "Point(7) = {XmaxCore, YmaxCore, ZminCore, mshC};\n"
out += "Point(8) = {XminCore, YmaxCore, ZminCore, mshC};\n"
out += "\n"
out += "Point(9) = {XminBB, YmaxBB, Zpos, mshBB };\n"
out += "Point(10) = {XminBB, YminBB, Zpos, mshBB };\n"
out += "Point(11) = {XmaxBB, YminBB, Zpos, mshBB };\n"
out += "Point(12) = {XmaxBB, YmaxBB, Zpos, mshBB };\n"
out += "\n"
out += "Point(17) = {XminBB, YmaxBB, ZminBB, mshBB };\n"
out += "Point(18) = {XminBB, YminBB, ZminBB, mshBB };\n"
out += "Point(19) = {XmaxBB, YminBB, ZminBB, mshBB };\n"
out += "Point(20) = {XmaxBB, YmaxBB, ZminBB, mshBB };\n"
out += "\n"
out += "Line(1) = {10, 9};\n"
out += "Line(3) = {10, 11};\n"
out += "Line(4) = {11, 12};\n"
out += "Line(5) = {12, 9};\n"
out += "Line(6) = {4, 3};\n"
out += "Line(7) = {3, 2};\n"
out += "Line(8) = {2, 1};\n"
out += "Line(9) = {1, 4};\n"
out += "Line(10) = {4, 8};\n"
out += "Line(11) = {9, 17};\n"
out += "Line(12) = {12, 20};\n"
out += "Line(13) = {3, 7};\n"
out += "Line(14) = {2, 6};\n"
out += "Line(15) = {6, 5};\n"
out += "Line(16) = {1, 5};\n"
out += "Line(17) = {5, 8};\n"
out += "Line(18) = {17, 18};\n"
out += "Line(19) = {18, 10};\n"
out += "Line(20) = {18, 19};\n"
out += "Line(21) = {19, 11};\n"
out += "Line(22) = {20, 17};\n"
out += "Line(23) = {20, 19};\n"
out += "Line(24) = {8, 7};\n"
out += "Line(25) = {7, 6};\n"
out += "Line Loop(1) = {6, 13, -24, -10};\n"
out += "Plane Surface(1) = {1};\n"
out += "Line Loop(2) = {10, -17, -16, 9};\n"
out += "Plane Surface(2) = {2};\n"
out += "Line Loop(3) = {15, 17, 24, 25};\n"
out += "Plane Surface(3) = {3};\n"
out += "Line Loop(4) = {13, 25, -14, -7};\n"
out += "Plane Surface(4) = {4};\n"
out += "Line Loop(5) = {14, 15, -16, -8};\n"
out += "Plane Surface(5) = {5};\n"
out += "Line Loop(6) = {19, 1, 11, 18};\n"
out += "Plane Surface(6) = {6};\n"
out += "Line Loop(7) = {22, 18, 20, -23};\n"
out += "Plane Surface(7) = {-7};\n"
out += "Line Loop(8) = {19, 3, -21, -20};\n"
out += "Plane Surface(8) = {-8};\n"
out += "Line Loop(9) = {23, 21, 4, 12};\n"
out += "Plane Surface(9) = {-9};\n"
out += "Line Loop(10) = {12, 22, -11, -5};\n"
out += "Plane Surface(10) = {10};\n"
out += "Line Loop(11) = {1, -5, -4, -3};\n"
out += "Line Loop(12) = {9, 6, 7, 8};\n"
out += "Plane Surface(11) = {12};\n"
out += "Plane Surface(12) = {11, 12};\n"

out += "// electrodes  (in the core)\n"
out += "k=newp;\n"
for s in range(NumElectrodes):
    #if electrode_in_use[electrodes_id[s]]:
        out += "Point(k+%s)={ %s, %s, %s, mshE};\n" % (
        s + 1, positions[s][0], positions[s][1], positions[s][2])
        out += "Point{k+%s} In Surface{11};\n" % (s + 1)
    #out += 'Physical Point("s%s")  = { k+%s } ;\n' % (s, i + 1)
out += 'Physical Surface("faces") = { 6,7,8,9,10};\n'
out += "Surface Loop(1) = {1,2,3,4,5,11};\n"
out += "Volume(1) = {-1};\n"
out += "Surface Loop(2) = {6, 8, 9, 7, 10, 2, 1, 12, 5, 4, 3};\n"
out += "Volume(2) = {2};\n"
out += 'Physical Volume("padding")  = { 2 } ;\n'
out += 'Physical Volume("core")  = { 1 } ;\n'
# out += "Surface Loop(2) = {6, -8, -9, -7, 10, 2, 1, 12, 5, 4, 3};\n"

GEOFN2 = args.geofile + ".geo"
MSHN3 = args.geofile + ".msh"

open(GEOFN2, 'w').write(out)
print("geometry has been written to %s" % GEOFN2)
rp = subprocess.run(["gmsh", "-3", "optimize_netgen", "-o", MSHN3, GEOFN2])
print(">> GMSH mesh file %s was generated." % MSHN3)

dts = []
dps = []
for i, s in enumerate(electrodes_id):
    if stationsFMT:
         dts.append(stationsFMT % s)
    else:
         dts.append(s)
    dps.append(positions[i])
domain = ReadGmsh(MSHN3, 3, diracPoints=dps, diracTags=dts, optimize=True)
domain.write(args.flyfile)
print("Mesh written to file %s" % (args.flyfile))
saveSilo(silofile, tags=makeTagMap(Function(domain)))
print("Mesh written to silo file %s.silo" % silofile)

# data records:
INPUT.readline()
NumRecords=int(INPUT.readline().split()[0])
print("Number of records found is ", NumRecords)
INPUT.readline()

DATAOUT=open(args.datafile, 'w')
NDATA=0
for i in range(NumRecords):
    s=INPUT.readline().split()
    A, B, M, N, V=int(s[0]), int(s[1]), int(s[2]), int(s[3]), float(s[4])
    if electrode_in_use[A] and electrode_in_use[B] and electrode_in_use[M] and electrode_in_use[N] :
        DATAOUT.write("%d, %d, %d, %d, %e\n" %( A, B, M, N, V) )
        NDATA+=1
DATAOUT.close()
print(f"Data with {NDATA} records was written to ", DATAOUT.name)

STATIONOUT=open(args.stationfile, 'w')
NDATA=0
for s in range(NumElectrodes):
        STATIONOUT.write("%d, %g, %g, %g\n" %( electrodes_id[s], positions[s][0], positions[s][1], positions[s][2]) )
        NDATA+=1
DATAOUT.close()
print(f"{NDATA} Stations were written to ", STATIONOUT.name)

# ... configuration file:
out = ""
for line in open(os.path.join(CDTEMPLATE, "config_template.txt"), 'r').readlines():
    if not line.startswith("#.") and not line.startswith("true_properties"):
        out += line

CHANGE = {
    "created": datetime.now().strftime("%d.%m.%Y %H:%M"),
    "project": os.path.splitext(os.path.basename(args.INP))[0],
    "meshfile": args.flyfile,
    "stationfile": args.stationfile,
    "stationsFMT": stationsFMT,
    "schedulefile": args.datafile,
    "datafile": args.datafile,
    "dipoleInjections": True,
    "dipoleMeasurements": True,
    "datacolumns": ['R' ],
    "alpha1": 1.,
    "core": ["core"],
    "faces": ["faces"]
}
FN = args.configfile + ".py"
open(FN, 'w').write(out.format(**CHANGE)+"true_properties=None\n")
print(f"Configuration file {FN} created.")
print(f"Please check the content of the file.")