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

parser = argparse.ArgumentParser(description='extracts a station from a gmsh geo and creates a fullwaver style schedule to create synthetic data.', epilog="version 1/2021")
parser.add_argument(dest='geofile', metavar='GEO', type=str, help='gmsh geo file with numElectrodesX, numElectrodesX and SpacingElectrodes specification. ')
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='python setting configuration')
parser.add_argument('--numDipoles', '-d', dest='numDipoles', metavar='numDipoles', type=int, default=999999999, help='number of charging dipole per extra N/S/E/W electrode.')

#parser.add_argument('--fullwaver', '-f', dest='fullwaver', action='store_true', default=False, help='create FullWaver data set')


args = parser.parse_args()
config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")


numElectrodesX, numElectrodesX, SpacingElectrodes = None, None, None
for line in open(args.geofile, 'r'):
    if line.strip().startswith("numElectrodesX"):
        numElectrodesX=int(line.split('=')[1][:-2])
    elif line.strip().startswith("numElectrodesY"):
        numElectrodesY=int(line.split('=')[1][:-2])
    elif line.strip().startswith("SpacingElectrodes"):
        SpacingElectrodes=float(line.split('=')[1][:-2])

print(f"Dimension of electrode array read from file {args.geofile}.")
print(f"numElectrodesX= {numElectrodesX}.")
print(f"numElectrodesY = {numElectrodesY}.")
print(f"SpacingElectrodes = {SpacingElectrodes}.")

if None in [numElectrodesX , numElectrodesY, SpacingElectrodes]:
    raise ValueError("Parameter missing.")

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
