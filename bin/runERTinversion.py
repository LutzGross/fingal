#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

sys.path.append(os.getcwd())
from fingal import ERTInversionH1, ERTInversionGauss, ERTInversionH2, ERTInversionGaussWithDiagonalHessian
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
from esys.escript.pdetools import MaskFromBoundaryTag, MaskFromTag
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
from esys.escript.minimizer import MinimizerLBFGS

import logging
from datetime import datetime
#setEscriptParamInt('TOO_MANY_LINES', 6000)

parser = argparse.ArgumentParser(description='driver to invert an ERT survey', epilog="version 16/11/2023")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=str, default="restart", help='reststart file name')
parser.add_argument('--restart', '-r',  dest='restart', action='store_true', default=False, help="start from restart file. RESTARTFN need to be set and exist.")
parser.add_argument('--savememory', '-s',  dest='savememory', action='store_true', default=False, help="try to save memory at the costs of CPU time.")
parser.add_argument('--nooptimize', '-n',  dest='nooptimize', action='store_true', default=False, help="Don't calibrated the value for config.sigma_ref before iteration starts.")
parser.add_argument('--test', '-t',  dest='testonly', action='store_true', default=False, help="stop after rescaling config.sigma_ref.")
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


logger.info("** This is an ERT inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
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
# set the reference conductivity:

# define region with fixed conductivity:
if config.fixed_region_tags and isinstance(config.fixed_region_tags , list):
    fixedm=MaskFromTag(domain, *tuple(config.fixed_region_tags))
    if len(config.fixed_region_tags)> 0:
        logger.info("Tags of regions with fixed property function : %s"%(config.fixed_region_tags) )
else:
    fixedm=None
# ... create mask where robin BC are applied or potential is set to zero, respectively.
if  config.use_robin_condition_in_model:
    mask_face = makeMaskForOuterSurface(domain, taglist=config.faces_tags)
else:
    mask_face = MaskFromBoundaryTag(domain, *config.faces_tags)

# create cost function:
logger.info(f"Regularization type = {config.regularization_order}.")
logger.info(f"Regularization_w1 = {config.regularization_w1}.")
if 'GAUSS' in config.regularization_order.upper():
    logger.info(f"Regularization_length_scale = {config.regularization_length_scale}.")

# initialize cost function:
if config.regularization_order == "H1":
    if config.use_robin_condition_in_model:
        assert m_ref is None, "ERTInversionH1WithRobinCondition does not support a reference property function yet."
        logger.info(f"Robin boundary conditions are applied in Forward model.")
        costf=ERTInversionH1WithRobinCondition(domain, data=survey,
                             sigma_0_ref=config.sigma0_ref, fixTop=config.fix_top,
                             w1=config.regularization_w1, useL1Norm=config.use_L1Norm, epsilonL1Norm=config.epsilon_L1Norm,
                             maskFixedProperty=fixedm, maskOuterFaces= mask_face, dataRTolDC= config.data_rtol,
                             pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                             useLogMisfitDC= config.use_log_misfit_DC, logger=logger.getChild("ERT-H1-Robin"))
    else:
        costf=ERTInversionH1(domain, data=survey,
                             sigma_0_ref=config.sigma0_ref, fixTop=config.fix_top,
                         w1=config.regularization_w1, useL1Norm=config.use_L1Norm, epsilonL1Norm=config.epsilon_L1Norm,
                         maskFixedProperty=fixedm, maskZeroPotential= mask_face, dataRTolDC= config.data_rtol, m_ref=m_ref,
                         pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                         useLogMisfitDC= config.use_log_misfit_DC, logger=logger.getChild("ERT-H1"))
    dM_init = Scalar(0.0, Solution(domain))

elif config.regularization_order == "H2":
    assert not config.use_robin_condition_in_model, "H2 Regularization does not support robin_condition in the model."
    costf=ERTInversionH2(domain, data=survey, save_memory = args.savememory,
                         sigma_0_ref=config.sigma0_ref, reg_tol=None, fixTop=config.fix_top,
                         w1=config.regularization_w1, maskZeroPotential= mask_face, dataRTolDC= config.data_rtol,  m_ref=m_ref,
                         pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                         useLogMisfitDC= config.use_log_misfit_DC, logger=logger.getChild("ERT-H2"))
    dM_init = Vector(0.0, Solution(domain))
elif config.regularization_order == "Gauss":
    assert not config.use_robin_condition_in_model, "Gauss Regularization does not support robin_condition in the model."
    costf = ERTInversionGauss(domain, data=survey,
                                    sigma_0_ref=config.sigma0_ref, fixTop=config.fix_top,
                                    w1=config.regularization_w1, length_scale = config.regularization_length_scale,
                                    maskZeroPotential=mask_face, dataRTolDC= config.data_rtol,
                                    pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,  m_ref=m_ref,
                                    logclip=config.clip_property_function,
                                    useLogMisfitDC=config.use_log_misfit_DC, logger=logger.getChild("ERT-Gauss"))
    dM_init = Data(0.0, (4,), Solution(domain))
elif config.regularization_order == "DGauss":
    assert not config.use_robin_condition_in_model, "DGauss Regularization does not support robin_condition in the model."
    costf = ERTInversionGaussWithDiagonalHessian(domain, data=survey,
                                                   sigma_0_ref=config.sigma0_ref, fixTop=config.fix_top,
                                                   w1=config.regularization_w1, length_scale = config.regularization_length_scale,
                                                   maskZeroPotential=mask_face, dataRTolDC= config.data_rtol,
                                                   pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                                                   logclip=config.clip_property_function,  m_ref=m_ref,
                                                   useLogMisfitDC=config.use_log_misfit_DC, logger=logger)
    dM_init = Data(0.0, (4,), Solution(domain))
else:
    raise ValueError("Unknown regularization type "+config.regularization_order)

# set up solver:
if not args.nooptimize:
    new_sigma_ref=costf.fitSigmaRef()
    logger.info(f"New value for config.sigma_ref = {new_sigma_ref}.")
    costf.updateSigma0Ref(new_sigma_ref)

if args.testonly:
    exit(0)
def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
    if args.RESTARTFN and iterCount >0:
        m.dump(args.RESTARTFN)
        logger.info(f"restart file {iterCount} for step {args.RESTARTFN} created.")
    #print(f"snapshot for step {k} saved.")
    #saveSilo("snapshot_"+args.OUTFILE, m=x)


# set up solver:
solver= MinimizerLBFGS(F=costf, iterMax=config.imax, logger=logger.getChild("LBFGS"))
solver.getLineSearch().setOptions(interpolationOrder=config.interpolation_order)
solver.setOptions(m_tol=config.m_tolerance, truncation=config.truncation, restart=config.restart, grad_tol=config.g_tolerance)
solver.setCallback(myCallback)
# initial solution
if args.restart:
    kk=[v for i,v in enumerate(os.listdir()) if v.startswith(args.RESTARTFN) ]
    if kk:
        dM_init=load(args.RESTARTFN, domain)
        logger.info(f"Restart file {args.RESTARTFN} read. Initial dM = {str(dM_init)}.")

        if config.regularization_order == 1:
            assert dM_init.getShape() == ()
        else:
            assert dM_init.getShape() == (4, )
# run solver:
solver.run(dM_init)
print(costf.getStatistics())
dM=solver.getResult()
m=costf.extractPropertyFunction(dM)
sigma=costf.getSigma0(m, applyInterploation=False)
if args.vtk:
    saveVTK(config.outfile, sigma=sigma, tag=makeTagMap(Function(domain)))
    logger.info(f"Result written to {config.outfile}.vtu.")
else:
    saveSilo(config.outfile, sigma=sigma, tag=makeTagMap(Function(domain)))
    logger.info(f"Result written to {config.outfile}.silo.")

if args.xyz:
    sigmas=interpolate(sigma, ReducedFunction(domain))
    X=sigmas.getFunctionSpace().getX()
    if isinstance(config.core, list):
        saveDataCSV(config.outfile+".csv", pos0=X[0], pos1=X[1], pos2=X[2], sigma =sigmas, mask=m)
        logger.info(f"Result written to {config.outfile}.csv. Core region is {config.core}.")
        del X, sigmas, m
    else:
        saveDataCSV(config.outfile+".csv", pos0=X[0], pos1=X[1], pos2=X[2], sigma =sigmas)
        logger.info(f"Tesult written to {config.outfile}.csv.")
        del X, sigmas
logger.info("All done - Have a nice day!")