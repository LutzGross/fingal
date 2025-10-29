#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
sys.path.append(os.getcwd())

import argparse
import numpy as np
from fingal import readElectrodeLocations, readSurveyData, ERTSensitivity
from esys.finley import ReadMesh, ReadGmsh
from esys.escript.pdetools import MaskFromBoundaryTag, MaskFromTag
from esys.weipa import saveVTK, saveSilo
#from esys.escript.pdetools import Locator, MaskFromTag

parser = argparse.ArgumentParser(description='creates ERT sensitivity map.')
parser.add_argument('--threshold', '-t', dest='threshold',  type=float, default=0.1, help='sensitivity threshold')
parser.add_argument('--sigma_true', '-s', dest='truesigma',  action='store_true', default=False, help='use the true conductivity distribution set in configuation if available.')
parser.add_argument('--vtk', '-v', dest='usevtk',  action='store_true', default=False, help='use the vtk output format.')
parser.add_argument('--obs', '-o',  dest='obs', metavar='OBS', type=str, default = None, help="request sensitivity for specific observation, format = `A,B,M,N`")
parser.add_argument('--file', '-f',  dest='outfile', metavar='OUTFILE', type=str, default='sensitivity', help="file for saving sensitivity for visualization (no extension) (default: `sensitivity`).")
parser.add_argument('--mesh', '-m', dest='msh', metavar='meshfile', type=str, default=None, help='name of mesh file (fly or msh). If not set the meshfile in the configuration file is used.')
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')

args = parser.parse_args()

print("** This creates a ERT sensitivity map **")
config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

if args.msh:
    if os.path.splitext(args.msh)[1] == ".msh":
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
    else:
        domain=ReadMesh(args.msh)
        print("Mesh read from file %s"%args.msh)
else:
    domain = ReadMesh(config.meshfile)
    print("Mesh read from file %s" % config.meshfile)

if config.schedulefile and os.path.exists(config.schedulefile):
    schedulefile = config.schedulefile
else:
    schedulefile = config.datafile
schedule = readSurveyData(schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates,
                        columns=[],
                        dipoleInjections=config.dipoleInjections,
                        dipoleMeasurements=config.dipoleMeasurements,
                        delimiter=config.datadelimiter,
                        commend='#', printInfo=True)
print(f"survey schedule read from {schedulefile}.")

out= {"tag" : makeTagMap(Function(domain)) }

maskZeroPotential = MaskFromBoundaryTag(domain, *config.faces_tags)
# -----------------------------------------------------------------------------------
runner=ERTSensitivity(domain, schedule,  sigma_src=config.sigma0_ref,
                    maskZeroPotential = maskZeroPotential,
                    stationsFMT=config.stationsFMT, printInfo = True)

if args.truesigma and callable(config.true_properties):
    sigma, M_n = config.true_properties(domain)
    print("True conductivity set in configuration file is used.")
    runner.setProperties(sigma=interpolate(sigma_0, Function(domain)))
    sigma_0.expand()
    out["sigma0"] =  sigma
else:
    print("Constant conductivity is used.")

if args.obs:
    obsS=runner.getSensitivityDensity(*tuple([int(ST) for ST in args.obs.split(',')]))
    kw="S"
    for ST in args.obs.split(','):
        kw+=f"_{int(ST)}"
    out[kw] = abs(obsS)

sensitivity=runner.getTotalSensitivityDensity()
resolution = (args.threshold / (sensitivity * runner.GEOMETRY_FACTOR) ) ** (1./3)
resolutionloss = (Lsup(sensitivity)/sensitivity)**(1./3)
out['Sensitivity'] = sensitivity
out['ResolutionLoss'] = resolutionloss
out['Resolution'] = resolution

if args.usevtk:
    saveVTK(args.outfile , **out)
    print(f'Values written to file {args.outfile}.vtk.')
else:
    saveSilo(args.outfile , **out)
    print(f'Values written to file {args.outfile}.silo.')
print(f"All done! Have a nice day.")