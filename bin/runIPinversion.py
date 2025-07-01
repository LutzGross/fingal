#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys
sys.path.append(os.getcwd())
from fingal import IPInversionH1, IPInversionH2, readElectrodeLocations, readSurveyData
from esys.finley import ReadMesh
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
from esys.escript.pdetools import MaskFromTag, MaskFromBoundaryTag
from esys.escript.minimizer import MinimizerLBFGS
#from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag



import logging
from datetime import datetime
#setEscriptParamInt('TOO_MANY_LINES', 6000)

parser = argparse.ArgumentParser(description='driver to invert an ERT survey', epilog="version 16/3/2018")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=None, default="restart", help='reststart file name')
parser.add_argument('--restart', '-r',  dest='restart', action='store_true', default=False, help="start from restart file. RESTARTFN need to be set and exist.")
parser.add_argument('--savememory', '-s',  dest='savememory', action='store_true', default=False, help="try to save memory at the costs of CPU time.")
parser.add_argument('--nooptimize', '-n',  dest='nooptimize', action='store_true', default=False, help="Don't calibrated the value for config.sigma_ref before iteration starts.")
parser.add_argument('--test', '-t',  dest='testonly', action='store_true', default=False, help="Calculates a new value for sigma_ref and then stops.")
parser.add_argument('--vtk', '-v',  dest='vtk', action='store_true', default=False, help="VTK format is used for output otherwise silo is used.")
parser.add_argument('--xyz', '-x',  dest='xyz', action='store_true', default=False, help="CSV file is create where results for the core region are written only.")
parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=False, help="shows more information.")
args = parser.parse_args()

logger=logging.getLogger('Inversion')
if args.debug:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)


config = importlib.import_module(args.config)


logger.info("** This is an IP inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
logger.info("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
logger.info("%s electrode locations read from %s."%(len(elocations), config.stationfile))

domain=ReadMesh(config.meshfile)
logger.info("Mesh read from "+config.meshfile)

survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns, 
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

# .... Load reference .... Work to do here!
m_ref = None
# set the reference conductivityand normalized chargeability:

# define region with fixed conductivity:
if config.fixed_region_tags and isinstance(config.fixed_region_tags , list):
    fixedm=MaskFromTag(domain, *tuple(config.fixed_region_tags))
    if len(config.fixed_region_tags)> 0:
        logger.info("Tags of regions with fixed property function : %s"%(config.fixed_region_tags) )
else:
    fixedm=None
# ... create mask where potential is set to zero:
mask_face = MaskFromBoundaryTag(domain, *config.faces_tags)
# create cost function:
logger.info(f"Regularization type = {config.regularization_order}.")
logger.info(f"Regularization_w1 = {config.regularization_w1}.")

# initialize cost function:
if config.regularization_order in ["H1_0", "H1"] :
    costf = IPInversionH1(domain, data=survey, sigma_0_ref=config.sigma0_ref, Mn_ref = config.Mn_ref,
                            w1 = config.regularization_w1, theta = config.regularization_theta,
                            maskZeroPotential=mask_face, maskFixedProperty=fixedm, fix_top=config.fix_top,
                            stationsFMT=config.stationsFMT, pde_tol= config.pde_tol,
                            weightingMisfitDC=config.regularization_weighting_DC_misfit,
                            useLogMisfitIP = config.use_log_misfit_IP, useLogMisfitDC=config.use_log_misfit_DC,
                            dataRTolDC = config.data_rtol, dataRTolIP = config.data_rtol,  m_ref=m_ref,
                            logclip=config.clip_property_function,
                            zero_mean_m=config.regularization_order == "H1_0",
                            logger = logger.getChild(f"IP-{config.regularization_order}"))
    dM_init = Data(0.0, (2,), Solution(domain))
elif config.regularization_order in ["H2_0", "H2"]:
    costf = IPInversionH2(domain, data=survey,  maskZeroPotential=mask_face,
                            pde_tol= config.pde_tol,
                            stationsFMT=config.stationsFMT, m_ref=m_ref,
                            useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC=config.data_rtol,
                            useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                            weightingMisfitDC=config.regularization_weighting_DC_misfit,
                            sigma_0_ref=config.sigma0_ref, Mn_ref= config.Mn_ref,
                            w1=config.regularization_w1, theta=config.regularization_theta,
                            fixTop=config.fix_top,  zero_mean_m = config.regularization_order == "H2_0",
                            logclip=config.clip_property_function, m_epsilon=1e-18,
                            length_scale=config.regularization_length_scale,
                            reg_tol=None, save_memory=args.savememory,
                            logger=logger.getChild(f"IP-{config.regularization_order}") )
    dM_init = Data(0.0, (6,), Solution(domain))
else:
    raise ValueError("Unknown regularization type " + config.regularization_order)

# set up solver:
if not args.nooptimize:
    new_sigma_ref, new_Mn_ref =costf.fitSigmaAndMnRef()
    logger.info(f"New value for config.sigma_ref = {new_sigma_ref}.")
    logger.info(f"New value for config.Mn_ref = {new_Mn_ref}.")
    costf.updateSigma0Ref(new_sigma_ref)
    costf.updateMnRef(new_Mn_ref)
if args.testonly:
    exit(0)

def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
    if args.RESTARTFN and iterCount >0:
        m.dump(args.RESTARTFN)
        logger.info("restart file %s for step %s created."%(iterCount, args.RESTARTFN))

# set up solver:
solver= MinimizerLBFGS(F=costf, iterMax=config.imax, logger=logger)
solver.getLineSearch().setOptions(interpolationOrder=config.interpolation_order)
solver.setOptions(m_tol=config.m_tolerance, truncation=config.truncation, restart=config.restart, grad_tol=config.g_tolerance)
solver.setCallback(myCallback)
# initial solution
if args.restart and args.restart:
   kk=[v for i,v in enumerate(os.listdir()) if v.startswith(args.RESTARTFN) ]
   if kk:
        dM_init=load(args.RESTARTFN, domain)
        txt=str(m_init)
        print("restart file %s read. initial M = %s"%(args.RESTARTFN, txt))

# run solver:
solver.run(dM_init)
m=costf.extractPropertyFunction(solver.getResult())

sigma0 = costf.getSigma0(m)
Mn = costf.getMn(m)

if args.vtk:
    saveVTK(config.outfile, sigma0=sigma0, Mn=Mn, tag=makeTagMap(Function(domain)))
    if getMPIRankWorld() == 0: 
        print("result written to %s.vtu"%config.outfile)
else:
    saveSilo(config.outfile, sigma0=sigma0, Mn=Mn, tag=makeTagMap(Function(domain)))
    if getMPIRankWorld() == 0: 
            print("result written to %s.silo"%config.outfile)

if args.xyz:
    sigma0s=interpolate(sigma0, ReducedFunction(domain))
    Mns = interpolate(Mn, sigma0s.getFunctionSpace())
    X=sigmas.getFunctionSpace().getX()
    if isinstance(config.core, list):
        m=insertTaggedValues(Scalar(0., sigma0s.getFunctionSpace()),  **{ t: 1 for t in config.core })
        saveDataCSV(config.outfile+".csv", d0=X[0], d1=X[1], d2=X[2], sigma0 =sigma0s, Mn =Mns, mask=m)
        if getMPIRankWorld() == 0: 
            print("result written to %s.csv. Core region is %s."%(config.outfile, config.core))
        del X, sigmas, m
    else:
        saveDataCSV(config.outfile+".csv", d0=X[0], d1=X[1], d2=X[2],  sigma0 =sigma0s, Mn =Mns)
        if getMPIRankWorld() == 0: 
            print("result written to %s.csv"%config.outfile)
        del X, sigmas
if getMPIRankWorld() == 0: 
    print("All done - Have a nice day!")
            
