#!/usr/bin/python3
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
from esys.escript import makeTagMap, ReducedFunction
from fingal import readElectrodeLocations


GEOFILE = "mine3D.geo"
MSHFILE = "mine.msh"
FLYFILE = 'mine.fly'
STATIONFILE = "stations.csv"

SILOFILE = "mesh"
SILOFILE = None

DELIMITER = ","

elocations=readElectrodeLocations(STATIONFILE, delimiter=DELIMITER)
print("%s electrode locations read from %s."%(len(elocations), STATIONFILE))


rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
# rp.check_returncode()
print(f"GMSH mesh file {MSHFILE} was generated.")


# convert into a FLY file:
domain = ReadGmsh(MSHFILE, 3,
                  diracPoints=[ elocations[s] for s in elocations ],
                  diracTags=[ f"s{s}" for s in elocations ],
                  optimize=True)
domain.write(FLYFILE)
print(f"mesh file {FLYFILE} was created.")

if SILOFILE:
    saveSilo(SILOFILE, tag=makeTagMap(ReducedFunction(domain)))
    print(f'Mesh written to file {SILOFILE}.silo.')