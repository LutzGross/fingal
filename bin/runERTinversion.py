#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys
sys.path.append(os.getcwd())
from fingal import PotentialERT, readElectrodeLocations, readSurveyData, makeTagMap
from esys.finley import ReadMesh, ReadGmsh
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
from esys.downunder import MinimizerLBFGS
from esys.escript.pdetools import MaskFromTag

#from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag



import logging
from datetime import datetime
#setEscriptParamInt('TOO_MANY_LINES', 6000)

parser = argparse.ArgumentParser(description='driver to invert an ERT survey', epilog="version 16/3/2018")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--optimize', '-o',  dest='optimize', action='store_true', default=False, help="Calibrated the value for sigma0 before iteration starts.")
parser.add_argument('--testSigma0', '-t',  dest='testSigma0', action='store_true', default=False, help="Calculates a new value for sigma0 and then stops.")
parser.add_argument('--addStations', '-a',  dest='addStations', action='store_true', default=False, help="station locations need to be added to mesh file as point elements.")
parser.add_argument('--query', '-q',  dest='query', type=str, default=None, help="give file name (if set) for file (A,B,M,N, observation, prediction) ")
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

if getMPIRankWorld() == 0:
    print("** This is an ERT  inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
    print("configuration "+args.config+" imported.")



elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
if getMPIRankWorld() == 0:
    print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

if config.stationsFMT:
     dts=[ config.stationsFMT%s for s in elocations]
else:
     dts=[ s for s in elocations]

if os.path.splitext(config.meshfile)[1] == ".msh" :
    domain=ReadGmsh(config.meshfile, 3, 
                    diracPoints=[ elocations[s] for s in elocations], diracTags=dts, optimize=True )
else:
    if args.addStations:
        domain=ReadMesh(config.meshfile, 
                        diracPoints=[ elocations[s] for s in elocations], diracTags=dts, optimize=True )
    else:
        domain=ReadMesh(config.meshfile)
if getMPIRankWorld() == 0:
    print("Mesh read from "+config.meshfile)

      
survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns, 
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

# set the reference conductivity:
if isinstance(config.sigma0, dict):
    sigma0=Scalar(config.sigma_background,Function(domain))
    for k in config.sigma0:
        sigma0.setTaggedValue(k, config.sigma0[k])
else:
    sigma0=config.sigma0
if getMPIRankWorld() == 0: 
    print("Reference conductivity sigma0 = %s"%(str(sigma0)))

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
        
# create cost function:
costf=PotentialERT(domain, data=survey, w0=config.w0, w1=config.w1, sigma0=sigma0, region_fixed=fixedm, stationsFMT=config.stationsFMT, alpha0=config.alpha0, alpha1=config.alpha1, weightLogDefect=config.weightLogDefect, logclip=config.clip_property_function)


if args.optimize or args.testSigma0:
    sigma_opt, f_opt,  defect = costf.optimizeSigma0()
    if getMPIRankWorld() == 0: 
         print("Better value for sigma0 = %s, correction factor %s, defect = %s"%(sigma_opt, f_opt,  defect))
         print("update your configuration file %s"%args.config)
         if args.testSigma0: 
            print("And goodbye")
         else:
            print("Sigma0 will be updated.") 
    if args.testSigma0: 
        sys.exit()
    else:
        costf.scaleSigma0(f_opt)

# test gradient:
if False:
    x=length(domain.getX())
    m=(x/Lsup(x))**3*4

    dm=(domain.getX()[0]+domain.getX()[1]+0.5*domain.getX()[2])/Lsup(x)/3*0.01

    J0=costf.getValue(m)
    print("J0=%e"%J0)	
    G=costf.getGradient(m)
    #print("gradient = %s"%str(G))
    Dex=costf.getDualProduct(dm,G)
    for k in range(4, 13):
        a=0.5**k
        J=costf.getValue(m+a*dm)
        D=(J-J0)/a
        print("XX %d, %e %e, %e %e, %e"%(k,J0, J, D, Dex, (D-Dex)/a) )

# set up solver:
solver=MinimizerLBFGS(J=costf, m_tol=config.tolerance, J_tol=config.tolerance*10, imax=config.imax)
solver.setOptions(interpolationOrder=config.interpolation_order, truncation=config.truncation, restart=config.restart)

# run solver:
m_init=Scalar(0., Solution(domain) ) 
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
            
