#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys
sys.path.append(os.getcwd())
from esys.finley import ReadMesh, ReadGmsh
from fingal import *

import argparse


parser = argparse.ArgumentParser(description='converts a gmsh msh file to finley fly file. Location of stations are merged.', epilog="l.gross@uq.edu.au, version 11/12/2020")
parser.add_argument('--silo', '-s', dest='silo',  type=str, default=None, help='name silo filet to show mesh while converting')
parser.add_argument(dest='msh', metavar='gmshFile', type=str, help='gmsh mesh file')
parser.add_argument(dest='config', metavar='configfile', type=str, help='configuration (no Ypadding extension)')
args = parser.parse_args()


print("** Gmsh to Finley converter **")

config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")


if hasattr(config, "electrodes") and config.electrodes:
    if config.stationsFMT:
        dts = [ config.stationsFMT%s for s in config.electrodes ]
    else:
        dts = [ s for s in config.electrodes]
    dps=[ (config.electrodes[s][0],config.electrodes[s][1],0.) for s in config.electrodes]
else:
    elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
    print("%s electrode locations read from %s."%(len(elocations), config.stationfile))
    dts = []
    dps = []
    for s in elocations:
        if config.stationsFMT:
            dts.append(config.stationsFMT%s)
        else:
            dts.append(s)
        dps.append(elocations[s])
domain=ReadGmsh(args.msh, 3, diracPoints=dps, diracTags=dts, optimize=True )
print("Gmsh msh file read from %s"%args.msh)

domain.write(config.meshfile)
print("Mesh written to file %s"%config.meshfile)

if args.silo:
    saveSilo(args.silo, tags=makeTagMap(Function(domain)))
    print("Mesh written to silo file %s.silo"%args.silo)
