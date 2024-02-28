#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
sys.path.append(os.getcwd())

import argparse
import numpy as np
from fingal import readElectrodeLocations, readSurveyData, IPSynthetic, makeMaskForOuterSurface
from esys.finley import ReadMesh
from esys.weipa import saveVTK, saveSilo
#from esys.escript.pdetools import Locator, MaskFromTag

parser = argparse.ArgumentParser(description='creates a synthetic ERT/IP survey data set', epilog="l.gross@uq.edu.au, version 9/8/2023")
parser.add_argument('--topdepth', '-t',  dest='topdepth', default=5, type=int, help="depth of top of the anomaly in term of mean electrode difference")
parser.add_argument('--radius', '-r',  dest='radius', default=6, type=int, help="radius of the anomaly in term of mean electrode difference")
parser.add_argument('--offset', '-o',  dest='offset', default=5, type=int, help="offset of the anomaly form the center in term of mean electrode difference")
parser.add_argument('--increase', '-i',  dest='increase', default=100, type=float, help="raise factor for sigma and Mn in anomaly")
parser.add_argument('--noise', '-n',  dest='noise', default=0., metavar='NOISE', type=int, help="%% of noise to be added. (default is 0) ")
parser.add_argument('--fullwaver', '-f', dest='fullwaver',  action='store_true', default=False, help='creates a fullwaver-style survey.')
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--silo', '-s',  dest='silofile', metavar='SILO', help="silo file for saving mesh file and property distributions for visualization (no extension) (output if set).")
args = parser.parse_args()

print("** This creates a synthetic survey data set from properties set by config.true_properties **")
config = importlib.import_module(args.config)

print("configuration "+args.config+" imported.")
print(f"Data columns to be generated: {config.datacolumns}")




domain=ReadMesh(config.meshfile)
print("mesh read from "+config.meshfile)

elocations = readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s." % (len(elocations), config.stationfile))

schedule = readSurveyData(config.schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates,
                        columns=[],
                        dipoleInjections=config.dipoleInjections,
                        dipoleMeasurements=config.dipoleMeasurements,
                        delimiter=config.datadelimiter,
                        commend='#', printInfo=True)

# -------------- simple setup of anomalies: -----------------------------------------
stationlocations = np.array( [ elocations[s].tolist() for s in elocations])
x_min, x_max=min(stationlocations[:,0]), max(stationlocations[:,0])
y_min, y_max=min(stationlocations[:,1]), max(stationlocations[:,1])
z_min, z_max=min(stationlocations[:,2]), max(stationlocations[:,2])

delectrodes = 1./(np.mean([ 1./np.linalg.norm(stationlocations[i]-stationlocations[i+1:], axis=1).min() for i in range(stationlocations.shape[0]-1)] ))
center = np.array( [ (x_min+x_max)/2, (y_min+y_max)/2, (z_min+z_max)/2 - (args.topdepth + args.radius) * delectrodes ] )
direction =  np.array( [ (x_min-x_max), (y_min-y_max), (z_min-z_max) ] )
direction*=-1/np.linalg.norm(direction)

print(f'x-range electrodes = {x_min} - {x_max}.')
print(f'y-range electrodes = {y_min} - {y_max}.')
print(f'z-range electrodes = {z_min} - {z_max}.')
print(f'center electrodes = {center}.')
print(f'direction = {direction}.')
print(f'mean electrode distance = {delectrodes}.')

core=insertTaggedValues(Scalar(0., ReducedFunction(domain)), **{ t: 1 for t in config.core })

sigma_0=Scalar(config.sigma_ref, core.getFunctionSpace())
X=sigma_0.getX()
anomaly_mask_0 = core * whereNonPositive(length( ( X - (center + direction * (args.offset + args.radius)  )) * [1,0,1]) - args.radius )
sigma_0+= (config.sigma_ref * args.increase - sigma_0 ) * anomaly_mask_0

M_n = Scalar(config.Mn_ref, core.getFunctionSpace())
anomaly_mask_Mn =  core * whereNonPositive(length( ( X - (center - direction * (args.offset + args.radius)  )) * [1,0,1]) - args.radius )
M_n+= (config.Mn_ref * args.increase - M_n ) * anomaly_mask_Mn

M_n_faces = config.Mn_ref
# -----------------------------------------------------------------------------------
runner=IPSynthetic(domain, schedule,  sigma_src=1,
                    mask_faces = makeMaskForOuterSurface(domain, taglist=config.faces),
                    stationsFMT=config.stationsFMT,
                    createSecondaryData=True,
                    createFieldData=args.fullwaver,  printInfo = True)


runner.setProperties(sigma_0=interpolate(sigma_0, Function(domain)), sigma_0_faces=config.sigma_ref,
                     M_n=interpolate(M_n, Function(domain)), M_n_faces=M_n_faces)


runner.write(config.datafile , datacolumns = config.datacolumns, addNoise = args.noise>0,
                            rel_error = args.noise, delimiter=config.datadelimiter,
                            usesStationCoordinates= config.usesStationCoordinates,
                            iFMT="%d", dFMT="%.5g", xFMT="%e")
# -----------------------------------------------------------------------------------
if args.silofile:
    kwargs = { "Mn" : M_n, "sigma0" : sigma_0, "tag" : makeTagMap(Function(domain)) }
    A=list(elocations.keys())[0]
    iA=schedule.getStationNumber(A)
    kwargs[f'VS_s{A}'] = runner.source_potential[iA]
    kwargs[f'dV0_s{A}'] = runner.potential_0[iA]
    kwargs[f'V0_s{A}'] = runner.potential_0[iA] + runner.source_potential[iA]
    kwargs[f'V2_s{A}'] = runner.potential_0[iA]-runner.potential_oo[iA]
    if runner.createFieldData:
        kwargs[f'ES_s{A}'] = runner.source_field[iA]
        kwargs[f'dE0_s{A}'] = runner.field_0[iA]
        kwargs[f'E0_s{A}'] = runner.field_0[iA] + runner.source_field[iA]
        kwargs[f'E2_s{A}'] = runner.field_0[iA] - runner.field_oo[iA]
    saveSilo(args.silofile , **kwargs)
    print(f'values {kwargs.keys()} written to file {args.silofile}.')
