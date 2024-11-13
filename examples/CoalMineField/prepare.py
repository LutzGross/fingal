#!/usr/bin/python3
from minetools import getGeometryFromGeoFile
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
from esys.escript import makeTagMap, ReducedFunction
import numpy as np

DATA_NORTH = '230911hf1.txt'
DATA_SOUTH = '230911ys1.txt'

GEOFILE = "mine.geo"
MSHFILE = "mine.msh"
SILOFILE = "mesh"

geo=getGeometryFromGeoFile(GEOFILE)
import config

data = {}
for name, fn, id_offset in [('south', DATA_SOUTH, 100), ('north', DATA_NORTH, 200)]:
    data_records=np.loadtxt(fn, skiprows=1, usecols = (0,1,2,3,6,7),
                dtype = np.dtype([('A',int), ('B', int), ('M', int), ('N',int), ("V", float), ('I', float)]) )

    A_min = min(min(data_records['A']), min(data_records['B']), min(data_records['N']), min(data_records['M']))
    A_max = max(max(data_records['A']), max(data_records['B']), max(data_records['N']), max(data_records['M']))
    n = 0
    for A, B, M, N, V, I in data_records:
        data[(id_offset+A-A_min+1,id_offset+B-A_min+1,id_offset+M-A_min+1,id_offset+N-A_min+1)] = V/I
        n+=1
    print(f"{name} - line:")
    print("\tnumber of records = ", n)
    print("\tnumber of electrodes = ", A_max-A_min+1)
    print("\tfirst electrode = ", A_min, '->', id_offset+1)
    print("\tlast electrode = ", A_max, '->', A_max + id_offset-A_min+1)


# extract some geometrical informationm from the geo file:

f=open(config.stationfile,'w')
for s in geo.Stations:
    f.write("%d, %s, %s, %s\n"%(s, geo.Stations[s][0], geo.Stations[s][1], geo.Stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(len(geo.Stations), config.stationfile ))

f=open(config.datafile,'w')
for A, B, M, N in data:
    f.write("%d, %d, %d, %d, %g\n"%(A, B, M, N, data[(A,B,M,N)]  ))
f.close()
print("%s observations  written to file %s"%(len(data), config.datafile ))

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