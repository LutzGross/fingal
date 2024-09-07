#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys
sys.path.append(os.getcwd())
from fingal import ERTInversion, readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
#from esys.downunder import MinimizerLBFGS
from esys.escript.pdetools import MaskFromTag
from esys.escript.minimizer import MinimizerLBFGS
#from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag



import logging
from datetime import datetime
#setEscriptParamInt('TOO_MANY_LINES', 6000)

parser = argparse.ArgumentParser(description='driver to invert an ERT survey', epilog="version 16/11/2023")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=None, default="restart", help='reststart file name')
parser.add_argument('--restart', '-r',  dest='restart', action='store_true', default=False, help="start from restart file. RESTARTFN need to be set and exist.")

parser.add_argument('--nooptimize', '-n',  dest='nooptimize', action='store_true', default=False, help="Don't calibrated the value for config.sigma_ref before iteration starts.")
parser.add_argument('--test', '-t',  dest='testonly', action='store_true', default=False, help="stop after rescaling config.sigma_ref.")
#parser.add_argument('--query', '-q',  dest='query', type=str, default=None, help="file name (if set) for output of (A,B,M,N, observation, prediction) ")
parser.add_argument('--vtk', '-v',  dest='vtk', action='store_true', default=False, help="VTK format is used for output otherwise silo is used.")
parser.add_argument('--xyz', '-x',  dest='xyz', action='store_true', default=False, help="CSV file is create where results for the core region are written only.")
parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=False, help="shows more information.")
args = parser.parse_args()

logger=logging.getLogger('inv')
if args.debug:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)


config = importlib.import_module(args.config)


print("** This is an ERT inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

domain=ReadMesh(config.meshfile)
print("Mesh read from "+config.meshfile)


survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns,
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

# set the reference conductivity:

# define region with fixed conductivity:
if config.fixed_region_tags and isinstance(config.fixed_region_tags , list):
    fixedm=MaskFromTag(domain, *tuple(config.fixed_region_tags))
    if len(config.fixed_region_tags)> 0:
        print("Tags of regions with fixed property function : %s"%(config.fixed_region_tags) )
else:
    fixedm=None
mask_face=makeMaskForOuterSurface(domain, taglist=config.faces_tags)

# create cost function:
costf=ERTInversion(domain, data=survey,
                   sigma_0_ref=config.sigma0_ref,
                   w1=config.w1, useL1Norm=config.use_L1Norm, epsilonL1Norm=config.epsilon_L1Norm,
                   mask_fixed_property=fixedm, mask_outer_faces = mask_face,
                   pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                   useLogMisfit= config.use_log_misfit_ERT, logger=logger)


# test gradient:
if True:
    #=====
    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    #pp=(x - inf(x))*(x - sup(x)) * (y - inf(y))* (y - sup(y))*(z - inf(z))
    pp = (z - inf(z))
    pp/=sup(abs(pp))
    #====

    x=length(domain.getX())
    m=x/Lsup(x)*pp
    ddm=(domain.getX()[0]+domain.getX()[1]+0.5*domain.getX()[2])/Lsup(x)/3*10*pp
    #ddm=pp*0.01
    print(str(m))
    print(str(ddm))
    dm=ddm
    m*=0
    args=costf.getArgumentsAndCount(m)
    G=costf.getGradientAndCount(m, *args)
    Dex = costf.getDualProductAndCount(dm, G)
    J0=costf.getValueAndCount(m,  *args)
    print("J0=%e"%J0)
    #print("gradient = %s"%str(G))

    print("XX log(a):\tJ0\t\tJ(a)\t\tnum. D\t\tD\t\terror O(a)\t\tO(1)")
    for k in range(4, 13):
        a=0.5**k
        J=costf.getValueAndCount(m+a*dm)
        D=(J-J0)/a
        print("XX \t%d:\t%e\t%e\t%e\t%e\t%e\t%e"%(k,J0, J, D, Dex, D-Dex, (D-Dex)/a) )
    1/0
# set up solver:
if not args.nooptimize:
    new_sigma_ref=costf.fitSigmaRef()
    print("new value for config.sigma_ref =",new_sigma_ref)
    costf.setSigma0Ref(new_sigma_ref)
    #costf.setSigmaSrc(new_sigma_ref)
if args.testonly:
    exit(0)
def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, norm_gradFm, args_m, failed):
    if args.RESTARTFN and iterCount >0:
        m.dump(args.RESTARTFN)
        if getMPIRankWorld() == 0:
            print("restart file %s for step %s created."%(iterCount, args.RESTARTFN))
    #print(f"snapshot for step {k} saved.")
    #saveSilo("snapshot_"+args.OUTFILE, m=x)


# set up solver:
solver= MinimizerLBFGS(F=costf, iterMax=config.imax, logger=logger)
solver.getLineSearch().setOptions(interpolationOrder=config.interpolation_order)
solver.setOptions(m_tol=config.m_tolerance, truncation=config.truncation, restart=config.restart, grad_tol=config.g_tolerance)
solver.setCallback(myCallback)
# initial solution
if args.restart:
   kk=[v for i,v in enumerate(os.listdir()) if v.startswith(args.RESTARTFN) ]
   if kk:
     m_init=load(args.RESTARTFN, domain)
     txt=str(m_init)
     if getMPIRankWorld() == 0:
             print("restart file %s read. initial M = %s"%(args.RESTARTFN, txt))
else:
    m_init = Scalar(0.0, Solution(domain))
# run solver:
solver.run(m_init)
m=solver.getResult()



sigma=costf.getSigma0(m)
if args.vtk:
    saveVTK(config.outfile, sigma=sigma, tag=makeTagMap(Function(domain)))
    if getMPIRankWorld() == 0:
        print("result written to %s.vtu"%config.outfile)
else:
    saveSilo(config.outfile, sigma=sigma, tag=makeTagMap(Function(domain)))
    if getMPIRankWorld() == 0:
            print("result written to %s.silo"%config.outfile)

if args.xyz:
    sigmas=interpolate(sigma, ReducedFunction(domain))
    X=sigmas.getFunctionSpace().getX()
    if isinstance(config.core, list):
        saveDataCSV(config.outfile+".csv", d0=X[0], d1=X[1], d2=X[2], sigma =sigmas, mask=m)
        if getMPIRankWorld() == 0:
            print("result written to %s.csv. Core region is %s."%(config.outfile, config.core))
        del X, sigmas, m
    else:
        saveDataCSV(config.outfile+".csv", d0=X[0], d1=X[1], d2=X[2], sigma =sigmas)
        if getMPIRankWorld() == 0:
            print("result written to %s.csv"%config.outfile)
        del X, sigmas
if getMPIRankWorld() == 0:
    print("All done - Have a nice day!")