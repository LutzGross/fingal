#!/usr/bin/python3
from minetools import getGeometryFromGeoFile
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
from esys.escript import makeTagMap, ReducedFunction
GEOFILE = "mine.geo"
MSHFILE = "mine.msh"
SILOFILE = "mesh"

import config
# extract some geometrical informationm from the geo file:
geo=getGeometryFromGeoFile(GEOFILE)
f=open(config.stationfile,'w')
for s in geo.Stations:
    f.write("%d, %s, %s, %s\n"%(s, geo.Stations[s][0], geo.Stations[s][1], geo.Stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(len(geo.Stations), config.stationfile ))

f=open(config.schedulefile,'w')
for A, B, M, N in geo.Schedule:
    f.write("%d, %d, %d, %d\n"%(A, B, M, N ))
f.close()
print("%s observations  written to file %s"%(len(geo.Schedule), config.schedulefile ))

# generate mesh:

rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
# rp.check_returncode()
print(f"GMSH mesh file {MSHFILE} was generated.")

# convert into a FLY file:
domain = ReadGmsh(MSHFILE, 3,
                  diracPoints=list(geo.Stations.values()),
                  diracTags=[ f"s{i}" for i in geo.Stations ],
                  optimize=True)
domain.write(config.meshfile)
print(f"mesh file {config.meshfile} was created.")

if SILOFILE:
    saveSilo(SILOFILE, tag=makeTagMap(ReducedFunction(domain)))
    print(f'Mesh written to file {SILOFILE}.silo.')