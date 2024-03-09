#!/usr/bin/python3
"""
script to set up a station and schedule file for parallel Wenner lines.
"""
import numpy as np
import argparse
import importlib, sys, os
from datetime import datetime
from fingal import makeWennerArray

parser = argparse.ArgumentParser(description='Generates a schedule and station file for parallel wenner lines along the x0 axis.', epilog="version 3/2024")
parser.add_argument('--numElectrodesPerLine', '-n', dest='numElectrodesPerLine', metavar='numElectrodesPerLine', type=int, default=32, help='number of electrodes per line (32).')
parser.add_argument('--numLines', '-l', dest='numLines', metavar='numLines', type=int, default=1, help='number of parallel lines (1)')
parser.add_argument('--lineLength', '-L', dest='lineLength', metavar='lineLength', type=float, default=320, help='length of lines (320)')
parser.add_argument('--lineDistance', '-D', dest='lineDistance', metavar='lineDistance', type=float, default=10, help='distance of of lines (10).')
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='python setting configuration')
args = parser.parse_args()

config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")


SurveyLengthX=args.lineLength
SurveyLengthY=args.lineDistance * ( args.numLines - 1)
#
k1=1;
stations={}
for i in range(args.numLines):
    for j in range(args.numElectrodesPerLine):
        xpos=(j/(args.numElectrodesPerLine-1) -0.5) * SurveyLengthX
        if args.numLines >1:
            ypos = (i/(args.numLines-1) -0.5) * SurveyLengthY
        else:
            ypos = 0.
        stations[k1+args.numElectrodesPerLine * i + j]=(xpos, ypos, 0.0)
assert len(stations) == args.numElectrodesPerLine * args.numLines, "generated and expected number of stations don't match"
print(f"{args.numLines} lines with {args.numElectrodesPerLine} generated.")
print(f"Length x width = {SurveyLengthX} x {SurveyLengthY}.")
print(f"electrode and line spacing = {args.lineLength/(args.numElectrodesPerLine-1)}, {args.lineDistance}.")

with open(config.stationfile, 'w') as fout:
    for id, X in stations.items():
        fout.write("%d%s %g%s %g%s %g\n" % (id, config.stationdelimiter, X[0], config.stationdelimiter, X[1], config.stationdelimiter, X[2]))
print(f'Station were written to {config.stationfile}')



cc=0
f=open(config.schedulefile,'w')
for i in range(args.numLines):
    k0=k1 + args.numElectrodesPerLine * i
    schedule=makeWennerArray(numElectrodes=args.numElectrodesPerLine, id0=0)
    for A,B,M, N in schedule:
        f.write("%s, %s, %s, %s\n"%(A+k0, B+k0,M+k0, N+k0))
        cc+=1
f.close()
print("%s observations (dipole-dipole Wenner) written to file %s"%(cc,f.name))

1/0
#-------------
stationsFMT="s%s"
locfile="stations.loc"
schedulefile="schedule.csv"
datafile="data.csv"
meshfile="dump.fly"
gmshfile="dump.msh"
silofile="mesh"

ID0=1
CDTEMPLATE="../../../bin"
TruePropertyDef="""
def true_properties(domain): 
    sigma_dump=1e-4     # for DC 
    sigma_target=1e-3   # for DC
    Mn_target=0.6/(1-0.6)*sigma_target

    from esys.escript import Scalar, Function, ReducedFunction, interpolate, whereNonPositive
    ThicknessHallow = {ThicknessHallow}
    TargetDepthFromCrown  = {TargetDepthFromCrown}
    TargetOffsetFromTop = {TargetOffsetFromTop}
    DumpHeightAboveDumpFoot = {DumpHeightAboveDumpFoot}
    DumpCrone = {DumpCrone}
    TopDumpEdgeFromFoot = {TopDumpEdgeFromFoot}

    assert ThicknessHallow <  TargetDepthFromCrown, "thickness of hallow layer needs to be less then the depth of target"
    X=ReducedFunction(domain).getX()
    mdump=Scalar(0., ReducedFunction(domain))
    mdump.setTaggedValue('Dump', 1)

    H=(DumpHeightAboveDumpFoot-TargetDepthFromCrown)/(DumpCrone + TopDumpEdgeFromFoot-TargetOffsetFromTop)*(X[0]-TargetOffsetFromTop) + TargetDepthFromCrown
    m1= whereNonPositive(X[2]+TargetDepthFromCrown)*whereNonPositive(X[2]+H)*mdump
    m1h =whereNonPositive(X[2] + TargetDepthFromCrown-ThicknessHallow)*whereNonPositive(X[2]+H-ThicknessHallow)*mdump

    sigma_true=Scalar(sigma_background, ReducedFunction(domain))
    sigma_true.setTaggedValue('Dump', sigma_dump)
    sigma_true=m1h*sigma_target+(1-m1h)*sigma_true

    Mn_true=Scalar(0., ReducedFunction(domain))
    #Mn_true.setTaggedValue('Dump', sigma_dump)
    Mn_true=m1*Mn_target+(1-m1)*Mn_true

    return interpolate(sigma_true, Function(domain)), interpolate(Mn_true, Function(domain))
"""
# === serach for information in geo file:
NumElectrodes, SurveyLineOffsetFromTopEdge, SurveyLineOffsetFromFootEdge, DumpAngleFoot, TargetOffsetFromTop, TargetDepthFromCrown, DumpHeightAboveDumpFoot, ThicknessHallow, DumpCrone = None, None, None, None, None, None, None, None, None

for line in open(args.geofile, 'r'):
    if line.find(';') > -1:
        lin2=line[0:line.index(';')]
        if lin2.strip().startswith("NumElectrodes"):
            NumElectrodes=int(lin2.split('=')[1])
        elif lin2.strip().startswith("SurveyLineOffsetFromTopEdge"):
            SurveyLineOffsetFromTopEdge=float(lin2.split('=')[1])
        elif lin2.strip().startswith("SurveyLineOffsetFromFootEdge"):
            SurveyLineOffsetFromFootEdge=float(lin2.split('=')[1])
        elif lin2.strip().startswith("DumpCrone"):
            DumpCrone=float(lin2.split('=')[1])
        elif lin2.strip().startswith("DumpAngleFoot"):
            DumpAngleFoot=float(lin2.split('=')[1])
        elif lin2.strip().startswith("DumpHeightAboveDumpFoot"):
            DumpHeightAboveDumpFoot=float(lin2.split('=')[1])
        elif lin2.strip().startswith("ThicknessHallow"):
            ThicknessHallow=float(lin2.split('=')[1])
        elif lin2.strip().startswith("TargetDepthFromCrown"):
            TargetDepthFromCrown=float(lin2.split('=')[1])
        elif lin2.strip().startswith("TargetOffsetFromTop"):
            TargetOffsetFromTop=float(lin2.split('=')[1])


print(f"Electrode are read from file {args.geofile}.")
print(f"numElectrodes= {NumElectrodes}.")
print(f"SurveyLineOffsetFromTopEdge = {SurveyLineOffsetFromTopEdge}.")
print(f"SurveyLineOffsetFromFootEdge = {SurveyLineOffsetFromFootEdge}.")
print(f"DumpCrone = {DumpCrone}.")
print(f"DumpAngleFoot = {DumpAngleFoot}.")
print(f"DumpHeightAboveDumpFoot = {DumpHeightAboveDumpFoot}.")
print(f"TargetDepthFromCrown = {TargetDepthFromCrown}.")
print(f"TargetOffsetFromTop = {TargetOffsetFromTop}.")


if None in [NumElectrodes, DumpAngleFoot, ThicknessHallow, TargetDepthFromCrown, TargetOffsetFromTop, DumpHeightAboveDumpFoot, SurveyLineOffsetFromTopEdge, SurveyLineOffsetFromFootEdge, DumpCrone]:
    raise ValueError("Parameter missing.")
ElectrodeSpacing=(DumpCrone-SurveyLineOffsetFromTopEdge-SurveyLineOffsetFromFootEdge)/(NumElectrodes-1)
print(f"SpacingElectrodes = {ElectrodeSpacing}.")
TopDumpEdgeFromFoot = DumpHeightAboveDumpFoot/np.tan(DumpAngleFoot/180*np.pi)

Stations={}
for num in range(NumElectrodes):
    Stations[ID0+num]=(SurveyLineOffsetFromTopEdge + num * ElectrodeSpacing, 0, 0)

# ..... create electrode location file ..... :
if args.mkLocations:
    with open(locfile, 'w') as fout:
        for id, X in Stations.items():
            fout.write("%d , %g, %g, %g\n"%(id, X[0], X[1], X[2]))
    print('station were written to ', locfile)

# ..... create schedule ..... :
if args.mkSchedule:

    if args.Wenner:
        schedule=makeWennerArray(numElectrodes=NumElectrodes, id0=ID0)
    else:
        schedule = makeZZArray(numElectrodes=NumElectrodes, id0=ID0)
    with open(schedulefile, 'w') as fout:
        for rec in schedule:
            fout.write("%d, %d, %d, %d\n"%tuple(rec))
    print("schedule with ",len(schedule)," records was written to ", schedulefile)




if args.mkConfig:
    out = ""
    for line in open(os.path.join(CDTEMPLATE, "config_template.txt"), 'r').readlines():
        if not line.startswith("#.") and not line.startswith("true_properties"):
            out += line

    CHANGE = {
        "created": datetime.now().strftime("%d.%m.%Y %H:%M"),
        "project": args.config,
        "meshfile": meshfile,
        "stationfile": locfile,
        "stationsFMT": stationsFMT,
        "schedulefile": schedulefile,
        "datafile": datafile,
        "dipoleInjections": True,
        "dipoleMeasurements": True,
        "datacolumns": ['R', 'ERR_R', 'ETA', 'ERR_ETA' ],
        "alpha1": 1.,
        "core" : ["Base", "Dump"],
        "faces" : ["OuterFaces"]
    }
    FORMATPROPDEF={
        "ThicknessHallow" : ThicknessHallow,
        "TargetDepthFromCrown" : TargetDepthFromCrown,
        "TargetOffsetFromTop" : TargetOffsetFromTop,
        "DumpHeightAboveDumpFoot" : DumpHeightAboveDumpFoot,
        "DumpCrone" : DumpCrone,
        "TopDumpEdgeFromFoot" : TopDumpEdgeFromFoot
    }
    FN = args.config + ".py"
    open(FN, 'w').write(out.format(**CHANGE)+TruePropertyDef.format(**FORMATPROPDEF))

    print(f"Configuration file {FN} created.")
    print(f"Please check the content of the file.")

if args.mkMesh:
    import subprocess
    from esys.escript import Function, makeTagMap
    from esys.finley import ReadGmsh
    from esys.weipa import saveSilo
    # now we are ready to generate the 3D mesh:
    rp=subprocess.run(["gmsh", "-3",  " -optimize_netgen", "-o", gmshfile, args.geofile])
    #print(rp.stderr)
    #rp.check_returncode()
    print(">> GMSH mesh file %s generated."%gmshfile)
    # ....
    dts = []
    dps = []
    for id, X in Stations.items():
         dts.append(stationsFMT % id)
         dps.append(X)
    domain = ReadGmsh(gmshfile, 3, diracPoints=list(Stations.values()), diracTags=[stationsFMT % s for s in  Stations.keys() ], optimize=True)
    print("Gmsh msh file read from %s" % gmshfile)
    domain.write(meshfile)
    print("Mesh written to file %s" % meshfile)
    saveSilo(silofile, tags=makeTagMap(Function(domain)))
    print("Mesh written to silo file %s.silo" % silofile)