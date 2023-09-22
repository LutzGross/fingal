#!/usr/bin/python3
"""
this is a simple script to extract the array electrode configuration from the gsmsh *.geo file and creates a 
a station file and survey schedule file as defined in the configuration file.

th script looks for numElectrodesX, numElectrodesX and SpacingElectrodes and adds a electrodes at east west north and south
that act as an electrodes for charging dipole in a fullwaver style survey. 
"""
import numpy as np
import argparse
import importlib, sys, os
from datetime import datetime
from fingal import makeZZArray, makeWennerArray

parser = argparse.ArgumentParser(description='extracts a station from a gmsh geo and creates a fullwaver style schedule to create synthetic data.', epilog="version 7/2023")
parser.add_argument('--geo', '-g', dest='geofile', metavar='GEO', default="DumpGeometry.geo", type=str, help='gmsh geo file with numElectrodes etc specification. ')
parser.add_argument('--mesh', '-m', dest='mkMesh', action='store_true', default=False,  help='Mesh file (fly) is generated.')
parser.add_argument('--loc', '-l', dest='mkLocations', action='store_true', default=False,  help='electrode location file is generated.')
parser.add_argument('--schedule', '-s', dest='mkSchedule', action='store_true', default=False,  help='schedule file is generated.')
parser.add_argument('--wenner', '-w', dest='Wenner', action='store_true', default=False,  help='Wenner survey is used.')
parser.add_argument('--config', '-c', dest='mkConfig', action='store_true', default=False,  help='make new config file')
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='python setting configuration')
#parser.add_argument('--numDipoles', '-d', dest='numDipoles', metavar='numDipoles', type=int, default=999999999, help='number of charging dipole per extra N/S/E/W electrode.')
args = parser.parse_args()

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