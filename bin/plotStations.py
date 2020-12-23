#!/usr/bin/python3
from __future__ import print_function
from esys.escript import *
import importlib, sys, os
from fingal import *
from esys.escript import unitsSI as U
import numpy as np
import argparse
import os.path
sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='creates plot of station/electrode locations', epilog="l.gross@uq.edu.au, version 18/4/2018")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--image', '-i',  dest='plotfile', type=str, default='locs.png', help="image of electrode locations.")
parser.add_argument('--nolabels', '-n',  dest='nolabels', action='store_true', default=False, help="show no station labels.")
parser.add_argument('--markersize', '-m',  dest='markersize', type=int, default=8, help="marker size.")
parser.add_argument('--fondsize', '-f',  dest='fsize', type=int, default=18, help="fond size.")
parser.add_argument('--scale', '-s',  dest='scale', type=int, default=0, help="axis scale")
parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=True, help="shows more information.")
args = parser.parse_args()

asx=10.**(-args.scale)


print("** This plots station locations **")

config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

fn=None
if os.path.isfile(config.schedulefile):
    fn=config.schedulefile
      
if os.path.isfile(config.datafile):
    fn=config.datafile
if fn is None:
    raise IOError("unable to find survey file.")   

survey=readGalvanicSurveyData(fn, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=[], dipoleInjections=config.dipoleInjections , dipoleMeasurements=config.dipoleMeasurements ,  delimiter=config.datadelimiter)
print("%s observations read from %s."%(survey.getNumObservations(), fn))

from matplotlib import pyplot 
fsize=args.fsize
pyplot.rc('font', size=fsize)
pyplot.rc('axes', labelsize=fsize)
pyplot.rc('axes', titlesize=fsize)
pyplot.rc('xtick', labelsize=fsize)
pyplot.rc('ytick', labelsize=fsize)

fig = pyplot.figure(figsize=(3.2, 3.2), dpi=400)
ax = fig.add_subplot(111)

ins=survey.getListOfInjectionStations()
obs=survey.getObservationElectrodes()

# distance between observation electrodes:
mindist=1e99
average=0
cc=0
for s in obs:
    x0=survey.getStationLocation(s)
    d= [ np.linalg.norm(survey.getStationLocation(t)-x0) for t in obs if not s == t ]
    mindist=min( mindist, min(d))
    average+=min(d)
    cc+=1

print("mean minimum electrode distance = %s"%(average/cc))
print("minimum electrode distance = %s"%mindist)

    
both=[ s for s in ins if s in obs]
ins2=[ s for s in ins if s not in both]
obs2=[ s for s in obs if s not in both]


xins=[ survey.getStationLocation(s)[0]*asx for s in ins2]
yins=[ survey.getStationLocation(s)[1]*asx for s in ins2]
xobs=[ survey.getStationLocation(s)[0]*asx for s in obs2]
yobs=[ survey.getStationLocation(s)[1]*asx for s in obs2]
xboth=[ survey.getStationLocation(s)[0]*asx for s in both]
yboth=[ survey.getStationLocation(s)[1]*asx for s in both]


pyplot.plot(xins,yins, 'bo', label='inject', markersize=args.markersize)
pyplot.plot(xobs,yobs, 'yo',  label='observations', markersize=args.markersize)
pyplot.plot(xboth,yboth, 'ko',  label='both', markersize=args.markersize)
if args.scale == 0:
    pyplot.xlabel("x [m]")
    pyplot.ylabel("y [m]")
elif args.scale == 3:
    pyplot.xlabel("x [km]")
    pyplot.ylabel("y [km]")
else:
    pyplot.xlabel("x [10^%s m]"%args.scale)
    pyplot.ylabel("y [10^%s m]"%args.scale)

if not args.nolabels:
   for s, x,y in zip(ins2, xins,yins):
       ax.annotate(str(s),xy=(x,y), fontsize=18)

   for s, x,y in zip(obs2, xobs,yobs):
       ax.annotate(str(s),xy=(x,y), fontsize=18)

   for s, x,y in zip(both, xboth,yboth):
       ax.annotate(str(s),xy=(x,y), fontsize=18)


print(len(obs2), " pure observation electrodes")
print(len(ins2), " pure injection electrodes")
print(len(both), "injection & observations electrodes")
pyplot.tight_layout()
  
pyplot.savefig(args.plotfile)
print("image saved to file ",args.plotfile)
