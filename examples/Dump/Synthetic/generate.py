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
parser.add_argument('--loc', '-l', dest='mkLocations', action='store_true', default=False,  help='electrode location file is to be generated.')
parser.add_argument('--schedule', '-s', dest='mkSchedule', action='store_true', default=False,  help='schedule file is to be generated.')
parser.add_argument('--wenner', '-w', dest='Wenner', action='store_true', default=False,  help='Wenner survey is is used.')
parser.add_argument('--config', '-c', dest='mkConfig', action='store_true', default=False,  help='make new config file')
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='python setting configuration')
parser.add_argument('--numDipoles', '-d', dest='numDipoles', metavar='numDipoles', type=int, default=999999999, help='number of charging dipole per extra N/S/E/W electrode.')
args = parser.parse_args()

locfile="stations.loc"
schedulefile="schedule.csv"
datafile="data.csv"
meshfile="dump.fly"
ID0=1001
CDTEMPLATE="../../../bin"
TruePrepertyDef="""
def true_properties(domain): 
    from esys.escript import Scalar, Function
    sigma_true=Scalar(sigma_background, Function(domain))
    Mn_true=Scalar(Mn_background, Function(domain))
    return sigma_true, Mn_true
"""
# === serach for information in geo file:
NumElectrodes, SurveyLineOffsetFromTopEdge, SurveyLineOffsetFromFootEdge, DumpCrone = None, None, None, None

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

print(f"Electrode are read from file {args.geofile}.")
print(f"numElectrodes= {NumElectrodes}.")
print(f"SurveyLineOffsetFromTopEdge = {SurveyLineOffsetFromTopEdge}.")
print(f"SurveyLineOffsetFromFootEdge = {SurveyLineOffsetFromFootEdge}.")
print(f"DumpCrone = {DumpCrone}.")

if None in [NumElectrodes, SurveyLineOffsetFromTopEdge, SurveyLineOffsetFromFootEdge, DumpCrone]:
    raise ValueError("Parameter missing.")
ElectrodeSpacing=(DumpCrone-SurveyLineOffsetFromTopEdge-SurveyLineOffsetFromFootEdge)/(NumElectrodes-1)
print(f"SpacingElectrodes = {ElectrodeSpacing}.")

# ..... create electrode location file ..... :
if args.mkLocations:
    with open(locfile, 'w') as fout:
        for num in range(NumElectrodes):
            fout.write("%d,%g, %g, %g\n"%(ID0+num, SurveyLineOffsetFromTopEdge + num * ElectrodeSpacing, 0, 0))
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
        "stationsFMT": "%s",
        "schedulefile": schedulefile,
        "datafile": datafile,
        "dipoleInjections": False,
        "dipoleMeasurements": False,
        "datacolumns": ['R', 'E' ],
        "alpha1": 1.,
        "core" : ["Base", "Dump"],
        "faces" : ["OuterFaces"]
    }
    FN = args.config + ".py"
    open(FN, 'w').write(out.format(**CHANGE)+TruePrepertyDef)
    print(f"Configuration file {FN} created.")
    print(f"Please check the content of the file.")

#====

1/0
# ==== create station location file
station_file=""



config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")



# expected number of electrodes: 

numElectrodes=numElectrodesX*numElectrodesY+4;
print(f"Total number of electrodes = {numElectrodes}.")


SurveyLengthX =(numElectrodesX-1)*SpacingElectrodes
SurveyLengthY =(numElectrodesY-1)*SpacingElectrodes
SurveyCenterX = 0.
SurveyCenterY = 0.
print(f"Area of the electrodes (without extra N/S/E/W electrodes): {SurveyLengthX} x {SurveyLengthY}.")


stations={}
k0=0;
stations[k0+1]=( ( (-1.)/(numElectrodesX-1)-0.5)*SurveyLengthX+SurveyCenterX, SurveyCenterY, 0.0);
stations[k0+2]=(  ( (numElectrodesX+0.)/(numElectrodesX-1)-0.5)*SurveyLengthX+SurveyCenterX, SurveyCenterY, 0.0);
stations[k0+3]=( SurveyCenterX, ((-1.)/(numElectrodesY-1) -0.5)*SurveyLengthY+SurveyCenterY, 0.0);
stations[k0+4]=( SurveyCenterX,  ((numElectrodesY+0.)/(numElectrodesY-1) -0.5)*SurveyLengthY+SurveyCenterY, 0.0);


k1=k0+5;
for i in range(numElectrodesY):
  for j in range(numElectrodesX):
        xpos=(j/(numElectrodesX-1) -0.5)*SurveyLengthX+SurveyCenterX;
        ypos=(i/(numElectrodesY-1) -0.5)*SurveyLengthY+SurveyCenterY;
        stations[k1+i*numElectrodesX+j]=(xpos, ypos, 0.0)
assert len(stations) == numElectrodes, "generated and expected number of stations don't match"

print("Measurements South-North = %s to %s m"%(min([stations[s][1] for s in stations if s>=k1]),max([stations[s][1] for s in stations if s>=k1])))
print("Measurements West-East = %s to %s m"%(min([stations[s][0] for s in stations if s>=k1]),max([stations[s][0] for s in stations if s>=k1])))

print("extra station west at  %s"%str(stations[k0+1]))
print("extra station east at  %s"%str(stations[k0+2]))
print("extra station north at  %s"%str(stations[k0+4]))
print("extra station south at  %s"%str(stations[k0+3]))

f=open(config.stationfile,'w')
for s in stations:
    f.write("%d, %s, %s, %s\n"%(s,stations[s][0], stations[s][1], stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(numElectrodes,config.stationfile ))

# new we create a schedule:
recorders=[]
for i in range(0,numElectrodesY):
    for j in range(0,numElectrodesX,1+i%2):
         M=k1+i*numElectrodesX+j
         recorders.append(M)
print("%s recording stations."%(len(recorders)))


charging=[]
for ia in range(1,numElectrodesY,2):
    for ja in range(1,numElectrodesX,2):
        A=k1+ia*numElectrodesX+ja
        charging.append(A)
print("%s potential charging stations."%(len(charging)))

schedule=[]
args.numDipoles=min(args.numDipoles, len(charging))
for B in [k0+1, k0+2, k0+3, k0+4]:
    As=np.random.permutation(charging)[:args.numDipoles]
    print("selected chargers for %s: %s"%(B, As))
    for A in As:
        for M in  recorders:
            schedule.append((A,B,M))
print("number of observations = %s"%len(schedule))


f=open(config.schedulefile,'w')
for A,B,M in schedule:
    f.write("%s, %s, %s\n"%(A, B,M))
f.close()
print("%s schedules observations written to file %s"%(len(schedule),config.schedulefile))
