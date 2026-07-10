#!/usr/bin/python3
import os, subprocess
from esys.finley import ReadGmsh

GEOFILE="./mine3D.geo"
MSHFILE="./mine.msh"
FLYFILE="./mine.fly"
dp = (44.1925003270153, 56.8285531244474, 1.5)

rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
print(">> GMSH mesh file %s was generated." % GEOFILE)


domain = ReadGmsh(MSHFILE, 3, diracPoints=[dp], diracTags=["S1"], optimize=True)
domain.write(FLYFILE)
print("Mesh written to file %s" % (FLYFILE) )