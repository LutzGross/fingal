#!/usr/bin/python3
from minetools import getGeometryFromGeoFile
import subprocess
from esys.finley import ReadGmsh
from esys.weipa import saveSilo
from esys.escript import makeTagMap, ReducedFunction
GEOFILE = "mine.geo"
MSHFILE = "mine.msh"
FLYFILE = "mine.fly"
SILOFILE = "mesh"

# extract some geometrical informationm from the geo file:
geo=getGeometryFromGeoFile(GEOFILE)

# generate mesh:

rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
# rp.check_returncode()
print(f"GMSH mesh file {MSHFILE} was generated.")

# convert into a FLY file:
domain = ReadGmsh(MSHFILE, 3,
                  diracPoints=list(geo.Stations.values()),
                  diracTags=[ f"s{i}" for i in geo.Stations ],
                  optimize=True)
domain.write(FLYFILE)
print(f"mesh file {MSHFILE} was created.")

saveSilo(SILOFILE, tag=makeTagMap(ReducedFunction(domain)))
print(f'Mesh written to file {SILOFILE}.silo.')