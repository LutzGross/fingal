#!/usr/bin/python3
"""
script to set up a station and schedule file for parallel Wenner lines.
"""
import numpy as np
import argparse
import importlib, sys, os
from datetime import datetime
from fingal import makeWennerArray

parser = argparse.ArgumentParser(description='Generates a schedule and station file for parallel wenner lines along the x0 axis.', epilog="version 3/2024")
parser.add_argument('--numElectrodesPerLine', '-n', dest='numElectrodesPerLine', metavar='numElectrodesPerLine', type=int, default=32, help='number of electrodes per line (32).')
parser.add_argument('--numLines', '-l', dest='numLines', metavar='numLines', type=int, default=1, help='number of parallel lines (1)')
parser.add_argument('--lineLength', '-L', dest='lineLength', metavar='lineLength', type=float, default=320, help='length of lines (320)')
parser.add_argument('--lineDistance', '-D', dest='lineDistance', metavar='lineDistance', type=float, default=10, help='distance of of lines (10).')
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='python setting configuration')
args = parser.parse_args()

config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")


SurveyLengthX=args.lineLength
SurveyLengthY=args.lineDistance * ( args.numLines - 1)
#
k1=1;
stations={}
for i in range(args.numLines):
    for j in range(args.numElectrodesPerLine):
        xpos=(j/(args.numElectrodesPerLine-1) -0.5) * SurveyLengthX
        if args.numLines >1:
            ypos = (i/(args.numLines-1) -0.5) * SurveyLengthY
        else:
            ypos = 0.
        stations[k1+args.numElectrodesPerLine * i + j]=(xpos, ypos, 0.0)
assert len(stations) == args.numElectrodesPerLine * args.numLines, "generated and expected number of stations don't match"
print(f"{args.numLines} lines with {args.numElectrodesPerLine} generated.")
print(f"Length x width = {SurveyLengthX} x {SurveyLengthY}.")
print(f"electrode and line spacing = {args.lineLength/(args.numElectrodesPerLine-1)}, {args.lineDistance}.")

with open(config.stationfile, 'w') as fout:
    for id, X in stations.items():
        fout.write("%d%s %g%s %g%s %g\n" % (id, config.stationdelimiter, X[0], config.stationdelimiter, X[1], config.stationdelimiter, X[2]))
print(f'Station were written to {config.stationfile}')



cc=0
f=open(config.schedulefile,'w')
for i in range(args.numLines):
    k0=k1 + args.numElectrodesPerLine * i
    schedule=makeWennerArray(numElectrodes=args.numElectrodesPerLine, id0=0)
    for A,B,M, N in schedule:
        f.write("%s, %s, %s, %s\n"%(A+k0, B+k0,M+k0, N+k0))
        cc+=1
f.close()
print("%s observations (dipole-dipole Wenner) written to file %s"%(cc,f.name))
print("all done.")

