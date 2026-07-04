#!/usr/bin/python3
"""
this script generates a series of meshes of the same geometry with a global refinement factor
the idea is to double the number of nodes between meshes.
"""
import os, subprocess
from esys.finley import ReadGmsh
meshdir = "meshes"
#for k in [1, 2, 4, 8, 16, 32, 64, 128] : #, 256, 512, 1024]:
for k in [ 256 ]:
    rf = k ** (-1/3.)
    geofile = f"dom{k}.geo"
    mshfile = f"tmp.msh"
    flyfile = f"dom{k}.fly"

    print(f"mesh {k}, refinement factor = {rf:g}, estm. Nodes = {10 * k}K ")

    with open('domain.tpl', 'r') as file:
        content = file.read()
    content = content.replace("%RFMFAC%", str(0.7 * rf) )
    with open(os.path.join(meshdir, geofile), 'w') as file:
        file.write(content)
    print(os.path.join(meshdir, geofile) + " generated.")

    rp = subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o",mshfile, os.path.join(meshdir, geofile)])
    print(">> GMSH mesh file %s was generated." % mshfile)

    domain = ReadGmsh(mshfile, 3, diracPoints=[ (0.,0.,0.)], diracTags=["S1"], optimize=True)
    domain.write(os.path.join(meshdir, flyfile))
    print("Mesh written to file %s" % (os.path.join(meshdir, flyfile)))