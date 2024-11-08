#!/usr/bin/python3
import subprocess
from fingal import makeWennerArray
import numpy as np
from esys.finley import ReadGmsh
GEOFILE = "brick.geo"
MSHFILE = "brick.msh"
FLYFILE = 'brick.fly'
STATIONSFILE = 'stations.csv'
SCHEDULEFILE =  'schedule.csv'

LineOffset = 0
NumElectrodes =  8
ElectrodeSpacing = 5
LineLength = (NumElectrodes-1) * ElectrodeSpacing
Stations={}
for k in range(NumElectrodes):
    stid = 100 + (k + 1)
    Stations[stid] = np.array([-LineLength/2 + k * ElectrodeSpacing, LineOffset, 0.])

f=open(STATIONSFILE,'w')
for s in Stations:
    f.write("%d, %s, %s, %s\n"%(s, Stations[s][0], Stations[s][1], Stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(len(Stations), STATIONSFILE))

Schedule = makeWennerArray(numElectrodes=NumElectrodes, id0=101)
f=open(SCHEDULEFILE,'w')
for A, B, M, N in Schedule:
    f.write("%d, %d, %d, %d\n"%(A, B, M, N ))
f.close()
print("%s observations  written to file %s"%(len(Schedule), SCHEDULEFILE ))

# generate mesh:
rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
# rp.check_returncode()
print(f"GMSH mesh file {MSHFILE} was generated.")

# convert into a FLY file:
domain = ReadGmsh(MSHFILE, 3,
                  diracPoints=list(Stations.values()),
                  diracTags=[ f"s{i}" for i in Stations ],
                  optimize=True)
domain.write(FLYFILE)
print(f"mesh file {FLYFILE} was created.")

