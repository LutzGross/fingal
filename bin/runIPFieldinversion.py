#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
from datetime import datetime
sys.path.append(os.getcwd())
from fingal import readElectrodeLocations, readSurveyData, makeTagField, DCInversionByFieldIntensity, getR2, ChargeabilityInversionByField
from esys.finley import ReadMesh
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
from esys.downunder import MinimizerLBFGS
from esys.escript.pdetools import MaskFromTag
import logging


parser = argparse.ArgumentParser(description='driver to invert an electric field intensity survey (aka Fullwaver). Measurements may not be dipoles.', epilog="l.gross@uq.edu.au, version Jan 2021")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration (no py extension)')
parser.add_argument('--sigmaOnly', '-s',  dest='sigmaonly', action='store_true', default=False, help="solve for conductivity only.")
parser.add_argument('--restart', '-r', dest = 'restart', action ='store_true', default=False, help='restart for chargeability from restart file with conductivity.')
parser.add_argument('--optimize', '-o',  dest='optimize', action='store_true', default=False, help="Calibrated the value for sigma0 before iteration starts.(ignored)")
parser.add_argument('--testSigma0', '-t',  dest='testSigma0', action='store_true', default=False, help="Calculates a new value for sigma0 and then stops.")
parser.add_argument('--vtk', '-v',  dest='vtk', action='store_true', default=False, help="VTK format is used for output otherwise silo is used.")
parser.add_argument('--xyz', '-x',  dest='xyz', action='store_true', default=False, help="XYZ file for conductivity is create.")
parser.add_argument('--useTrueSigma', '-u',  dest='truesigma', action='store_true', default=False, help="Use true sigma for chargeability inversion.")
parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=False, help="shows more information.")
args = parser.parse_args()

logger=logging.getLogger('inv')
if args.debug:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

config = importlib.import_module(args.config)
if getMPIRankWorld() == 0:
    print("** This is an ERT/IP inversion using field intensity data @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
    print("configuration "+args.config+" imported.")


domain=ReadMesh(config.meshfile)
if getMPIRankWorld() == 0:
    print("mesh read from "+config.meshfile)

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
if getMPIRankWorld() == 0:
    print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns, 
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=False, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

# define region with fixed conductivity:
if isinstance(config.region_fixed , list):
    fixedm=MaskFromTag(domain, *tuple(config.region_fixed))
    if len(config.region_fixed)> 0 and getMPIRankWorld() == 0: 
        print("Tags of regions with fixed property function : %s"%(config.region_fixed) )
else:
    x=domain.getX()
    fixedm=whereZero(x[2]-inf(x[2]))
    del x
    if getMPIRankWorld() == 0: 
        print("Properties are fixed at the bottom of the domain.")

# define true properties if available:
if not config.true_properties is None:
    sigma_true, gamma_true = config.true_properties(domain)
    txt1, txt2 = str(sigma_true).replace("\n", ";"), str(gamma_true).replace("\n", ";")
    if getMPIRankWorld() == 0:
        print("True properties set:\n\tsigma = %s\n\tgamma = %s"%(txt1, txt2))
else:
    assert args.truesigma, 'If true conductivity is used true_properties must be set.'


if not args.restart and not args.truesigma:
    if getMPIRankWorld() == 0:
        print("**** Start field inversion: w0=%s, w1=%s, alpha0=%s, alpha1=%s"%(config.w0 , config.w1 , config.alpha0, config.alpha1))

    # set the reference conductivity:
    if isinstance(config.sigma0, dict):
        sigma0=Scalar(config.sigma_background,Function(domain))
        for k in config.sigma0:
            sigma0.setTaggedValue(k, config.sigma0[k])
    else:
        sigma0=config.sigma0
    if getMPIRankWorld() == 0: 
        txt1=str(sigma0)
        print("Reference conductivity sigma0 = %s"%(txt1))
        
    # create an instance of the 
    costf=DCInversionByFieldIntensity(domain, data=survey, L_stations=config.L_stations, w0=config.w0, w1=config.w1, 
                                  alpha0=config.alpha0, alpha1=config.alpha1, 
                                  sigma0=sigma0, region_fixed=fixedm, 
                                  stationsFMT=config.stationsFMT, 
                                  adjustStationLocationsToElementCenter=config.adjustStationLocationsToElementCenter,
                                  weightLogDefect=config.weightLogDefect, logclip=config.clip_property_function)
    # this is just for testing:
    if False:
        x=length(domain.getX())
        m=x/Lsup(x)*0.01
        #m=Scalar(0.2, Solution(domain) ) 
        dm=-abs(domain.getX()[0]/Lsup(x))**2
        args=costf.getArguments(m)
        J0=costf.getValue(m,*args)
        print("J0=",J0)
        G=costf.getGradient(m,*args)


        for k in range(1, 15):
            a=0.5**k
            J1=costf.getValue(m+a*dm)
            D=(J1-J0)/a
            Dex=costf.getDualProduct(dm,G)
            print("XX", a, "->", J1,J0,  D, Dex, (D-Dex), '->', (D-Dex)/a,)
        1/0   
    
    
    # correct sigma0 if required:
    if args.optimize or args.testSigma0:
        sigma_opt, f_opt,  defect = costf.optimizeSigma0()
        if getMPIRankWorld() == 0: 
            print("Better value for sigma0 = %s, correction factor %s, defect = %s"%(sigma_opt, f_opt,  defect))
            print("update your configuration file %s"%args.config)
            if args.testSigma0: 
                print("And goodbye!")
            else:
                print("Sigma0 will be updated.") 
        if args.testSigma0: 
            sys.exit()
        else:
            costf.scaleSigma0(f_opt)
            
    # set up solver:
    solver=MinimizerLBFGS(J=costf, m_tol=config.tolerance, J_tol=config.tolerance*10, imax=config.imax)
    solver.setOptions(interpolationOrder=config.interpolation_order, truncation=config.truncation, restart=config.restart)

    # run solver:
    m_init=Scalar(0., Solution(domain) ) 
    solver.run(m_init)
    m=solver.getResult()            
    solver.logSummary()
    sigma=costf.getSigma(m)

    del costf, solver
 
    if not config.true_properties is None:
        sigmai=interpolate(sigma, Function(domain))   
        R2=getR2(log(sigma_true), log(sigmai))
        if getMPIRankWorld() == 0:
            print("R^2 for log(conductivity):")
            print("\ttotal domain = %s"%R2)
        if config.core:
            R2=getR2(log(sigma_true), log(sigmai), chi=insertTaggedValues(Scalar(0., sigmai.getFunctionSpace()), **{ t: 1 for t in config.core }))
            if getMPIRankWorld() == 0: print("\tcore = %s"%R2)


        
    if config.restartfile is not None:
        sigma.dump(config.restartfile)
        if getMPIRankWorld() == 0: print("restart file %s for conductivity created."%config.restartfile)
        
else:
    if args.restart:
        sigma=load(config.restartfile, domain)
        if getMPIRankWorld() == 0: print(">>restart from file %s."%config.restartfile)
    else:
        if not config.true_properties is None:
            sigma=sigma_true.copy()
            if getMPIRankWorld() == 0: print(">>true sigma is used.") 
        else:
            raise ValueError("No sigma defined.")
        
# now we start the inversion for chargeability
txt=str(sigma)
if getMPIRankWorld() == 0: print("sigma = ",txt)

if not args.sigmaonly:
    print("**** Start chargeability inversion: w0=%s, w1=%s, alpha0=%s, alpha1=%s"%(config.w0gamma, config.w1gamma, config.alpha0gamma, config.alpha1gamma ))
    costf=ChargeabilityInversionByField(domain,data=survey, L_stations=config.L_stations, w0=config.w0gamma, w1=config.w1gamma, 
                                                      alpha0=config.alpha0gamma, alpha1=config.alpha1gamma, gamma0=config.eta0/(1-config.eta0), 
                                                      sigma=sigma, region_fixed=fixedm, 
                                                      stationsFMT=config.stationsFMT, weightLogDefect=config.weightLogDefect, logclip=config.clip_property_function,
                                                      adjustStationLocationsToElementCenter=config.adjustStationLocationsToElementCenter)
    # this is just for testing:
    if False:
        x=length(domain.getX())
        m=x/Lsup(x)*0.01
        #m=Scalar(0.2, Solution(domain) ) 
        dm=-abs(domain.getX()[0]/Lsup(x))**2*100
        args=costf.getArguments(m)
        J0=costf.getValue(m, *args)
        print("J0=",J0)
        G=costf.getGradient(m,*args)


        for k in range(1, 15):
            a=0.5**k
            J1=costf.getValue(m+a*dm)
            D=(J1-J0)/a
            Dex=costf.getDualProduct(dm,G)
            print("XX", a, "->", J1,J0,  D, Dex, (D-Dex), '->', (D-Dex)/a,)
        1/0   
        

    # set up solver:
    solver=MinimizerLBFGS(J=costf, m_tol=config.tolerance, J_tol=config.tolerance*10, imax=config.imax)
    solver.setOptions(interpolationOrder=config.interpolation_order, truncation=config.truncation, restart=config.restart)

    # run solver:
    m_init=Scalar(0., Solution(domain) ) 
    solver.run(m_init)
    m=solver.getResult()            
    solver.logSummary()    
    
    gamma=costf.getGamma(m)
    eta=gamma/(1+gamma)
    txt1=str(gamma)
    txt2=str(eta)
    if getMPIRankWorld() == 0: print(f"gamma = {txt1}, eta={txt2}")

    if not config.true_properties is None:

        gammai=interpolate(gamma, ReducedFunction(domain))  

        R2=getR2(gamma_true, gammai)
        R2log=getR2(log(gamma_true), log(gammai))


        print("R^2 modified chargeability gamma:")        
        print("\ttotal domain  = %s"%R2)
        print("\tlog, total domain = %s"%R2log)

        if config.core:
            mask=insertTaggedValues(Scalar(0., gammai.getFunctionSpace()), **{ t: 1 for t in config.core })
            R2=getR2(gamma_true, gammai, chi=mask)
            R2log=getR2(log(gamma_true), log(gammai), chi=mask)
            if getMPIRankWorld() == 0: 
                print("\tcore = %s"%R2)
                print("\tlog, core = %s"%R2log)
                
                

# assemblage of output: 
sigma.expand()
outargs={ "tag" : makeTagField(Function(domain)), "sigma":sigma } 

if not config.true_properties is None:
    outargs['sigma_error']=sigma/sigma_true
    sigma_true.expand()
    outargs['sigma_true']=sigma_true        

if not args.sigmaonly:
    outargs['gamma']=gamma
    if not config.true_properties is None:
        outargs['gamma_error']=gamma/gamma_true
        gamma_true.expand()
        outargs['gamma_true']=gamma_true
        
if args.vtk:
    saveVTK(config.outfile, **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.vtu"%config.outfile)
else:
    saveSilo(config.outfile,  **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.silo"%config.outfile)

if args.xyz:
    sigmas=interpolate(sigma, ReducedFunction(domain))
    X=sigmas.getFunctionSpace().getX()
    args={ "d0": X[0], "d1": X[1], "d2" :X[2], "sigma" : sigmas}
    if not args.sigmaonly:
        args['gamma']=interpolate(gamma, sigmas.getFunctionSpace())
    saveDataCSV(config.outfile+".csv", **args)
    if getMPIRankWorld() == 0: print("result written to %s.csv"%config.outfile)

if getMPIRankWorld() == 0: 
    print("All done - Have a nice day!")
