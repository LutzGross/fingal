#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys
sys.path.append(os.getcwd())
from fingal import IPInversion, readElectrodeLocations, readSurveyData
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

parser = argparse.ArgumentParser(description='driver to invert an ERT survey', epilog="version 16/3/2018")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=None, default="restart", help='reststart file name')
parser.add_argument('--restart', '-r',  dest='restart', action='store_true', default=False, help="start from restart file. RESTARTFN need to be set and exist.")

parser.add_argument('--optimize', '-o',  dest='optimize', action='store_true', default=False, help="Calibrated the value for sigma_ref before iteration starts.")
parser.add_argument('--test', '-t',  dest='testsigma_ref', action='store_true', default=False, help="Calculates a new value for sigma_ref and then stops.")
parser.add_argument('--query', '-q',  dest='query', type=str, default=None, help="file name (if set) for output of (A,B,M,N, observation, prediction) ")
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


print("** This is an IP inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

domain=ReadMesh(config.meshfile)
print("Mesh read from "+config.meshfile)

      
survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns, 
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

# set the reference conductivity:
if isinstance(config.sigma_ref, dict):
    sigma_ref=Scalar(config.sigma_background,Function(domain))
    for k in config.sigma_ref:
        sigma_ref.setTaggedValue(k, config.sigma_ref[k])
else:
    sigma_ref=config.sigma_ref
print("Reference conductivity sigma_ref = %s"%(str(sigma_ref)))

if config.Mn_ref and isinstance(config.Mn_ref, dict):
    Mn_ref=Scalar(0.,Function(domain))
    for k in config.Mn_ref:
        Mn_ref.setTaggedValue(k, config.sigma_ref[k])
else:
    Mn_ref=config.Mn_ref
print("Reference normalized chargeabilty Mn_ref = %s"%(str(Mn_ref)))
print("Background conductivity sigma_background = %s"%(str(config.sigma_background)))


# define region with fixed conductivity:
if config.region_fixed and isinstance(config.region_fixed , list):
    fixedm=MaskFromTag(domain, *tuple(config.region_fixed))
    if len(config.region_fixed)> 0:
        print("Tags of regions with fixed property function : %s"%(config.region_fixed) )
else:
    x=domain.getX()
    fixedm=whereZero(x[2]-inf(x[2]))
    del x
    print("Properties are fixed at the bottom of the domain.")


# create cost function:
costf=IPInversion(domain, data=survey,
                  sigma_0_ref=sigma_ref, Mn_ref=Mn_ref, sigma_background=config.sigma_background,
                  w1=config.w1, theta=config.theta, usel1=config.usel1, epsl1=config.epsl1,
                  mask_fixed_property=fixedm,
                  weighting_misfit_DC=config.weighting_misfit_ERT, pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                  logclip=config.clip_property_function,
                  EPSILON=1e-8, logger=logger)


#if args.optimize or args.testsigma_ref:
#    sigma_opt, f_opt,  defect = costf.optimizesigma_ref()
#    if getMPIRankWorld() == 0:
#         print("Better value for sigma_ref = %s, correction factor %s, defect = %s"%(sigma_opt, f_opt,  defect))
#         print("update your configuration file %s"%args.config)
#         if args.testsigma_ref:
#            print("And goodbye")
#         else:
#            print("sigma_ref will be updated.")
#    if args.testsigma_ref:
#        sys.exit()
#    else:
#        costf.scalesigma_ref(f_opt)

# test gradient:
if False:
    x=length(domain.getX())
    m=(x/Lsup(x))*[1,1]
    m+=(domain.getX()[1]/Lsup(x))*[1,0]

    ddm=(domain.getX()[0]+domain.getX()[1]+0.5*domain.getX()[2])/Lsup(x)/3*10
    print(str(m))
    print(str(ddm))
    dm=ddm*[1,0]
    J0=costf.getValueAndCount(m)
    print("J0=%e"%J0)	
    G=costf.getGradientAndCount(m)
    #print("gradient = %s"%str(G))
    Dex=costf.getDualProductAndCount(dm,G)
    print("XX log(a):\tJ0\t\tJ(a)\t\tnum. D\t\tD\t\terror O(a)\t\tO(1)")
    for k in range(4, 13):
        a=0.5**k
        J=costf.getValueAndCount(m+a*dm)
        D=(J-J0)/a
        print("XX \t%d:\t%e\t%e\t%e\t%e\t%e\t%e"%(k,J0, J, D, Dex, D-Dex, (D-Dex)/a) )

# set up solver:

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
if args.restart and args.restart:
   kk=[v for i,v in enumerate(os.listdir()) if v.startswith(args.RESTARTFN) ]
   if kk:
     m_init=load(args.RESTARTFN, domain)
     txt=str(m_init)
     if getMPIRankWorld() == 0:
             print("restart file %s read. initial M = %s"%(args.RESTARTFN, txt))
else:
    m_init = Data(0.,(2,), Solution(domain))
# run solver:
solver.run(m_init)
m=solver.getResult()


if args.query:
    sigma, potentials, dV=costf.getArguments(m)
    if getMPIRankWorld() == 0: 
        f=open(args.query,'w')
        for t in survey.tokenIterator(): # t=(A,B,M,N) (or so)
            f.write("%s, %e, %e\n"%(str(t)[1:-1], survey.getResistenceData(t), dV[t]))
        f.close()
        print("Data query/result file %s written."%args.query)
else:
    sigma=costf.getSigma(m)

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
        m=insertTaggedValues(Scalar(0., sigmas.getFunctionSpace()), **{ t: 1 for t in config.core })
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
            
