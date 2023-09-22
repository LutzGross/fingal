#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
from esys.finley import ReadMesh, ReadGmsh
from esys.weipa import saveVTK, saveSilo
import argparse
from fingal import makeTagMap


parser = argparse.ArgumentParser(description='creates a silo/vtk file from fly.')
parser.add_argument(dest='meshfile', metavar='meshfile', type=str, help='fly file')
parser.add_argument('--silo', '-s',  dest='silo', metavar='SILO', help="silo file to write.")
parser.add_argument('--vtk', '-v',  dest='vtk', metavar='SILO', help="vtk file to write.")
args = parser.parse_args()



domain=ReadMesh(args.meshfile)
print("mesh read from "+args.meshfile)

if args.silo is not None:
    saveSilo(args.silo,tag=makeTagMap(ReducedFunction(domain)))
    print(args.silo+".silo with tags has been generate")
if args.vtk is not None:
    saveSilo(args.vtk,tag=makeTagMap(ReducedFunction(domain)))
    print(args.vtk+".vtk with tags has been generated.")
