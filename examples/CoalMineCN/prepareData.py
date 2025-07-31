#!/usr/bin/python3
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
from esys.escript import makeTagMap, ReducedFunction
from fingal import readElectrodeLocations
import numpy as np

DATA1 = 'Data/230911hf1.txt', 200
DATA2 = 'Data/230911ys1.txt', 100

GEOFILE = "mine3D.geo"
MSHFILE = "mine.msh"
SILOFILE = "mesh"

#geo=getGeometryFromGeoFile(GEOFILE)
import config

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

data = {}
for fn, id_offset in [DATA1, DATA2]:
    data_records=np.loadtxt(fn, skiprows=1, usecols = (0,1,2,3,6,7),
                dtype = np.dtype([('A',int), ('B', int), ('M', int), ('N',int), ("V", float), ('I', float)]) )

    n = 0
    for A, B, M, N, V, I in data_records:
        data[(id_offset+A,id_offset+B,id_offset+M,id_offset+N)] = V/I
        n+=1
    print(f"line: {fn}:")
    print("\tid offset = ", id_offset)
    print("\tnumber of records = ", n)



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
                  diracPoints=[ elocations[s] for s in elocations ],
                  diracTags=[ f"s{s}" for s in elocations ],
                  optimize=True)
domain.write(config.meshfile)
print(f"mesh file {config.meshfile} was created.")

if SILOFILE:
    saveSilo(SILOFILE, tag=makeTagMap(ReducedFunction(domain)))
    print(f'Mesh written to file {SILOFILE}.silo.')