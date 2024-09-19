#!/usr/bin/python3

import importlib, sys, os
from esys.escript import *


from datetime import datetime
from esys.escript.minimizer import MinimizerLBFGS
sys.path.append(os.getcwd())
from fingal import readElectrodeLocations, readSurveyData, InterpolatorWithExtension, mapToDomain2D, createBackgroundTemperature
from fingal import ConductivityModelByTemperature
from esys.finley import ReadMesh
import numpy as np
from esys.weipa import saveSilo
import argparse
from esys.escript.pdetools import MaskFromTag, MaskFromBoundaryTag
import logging
#from IPtemperatureinversion import ConductivityModelByTemperature, createBackgroundTemperature, TemperatureBasedIPInversion, InterpolatorWithExtension, mapToDomain2D


parser = argparse.ArgumentParser(description='driver to invert an electric field intensity survey (aka Fullwaver). Measurements may not be dipoles.', epilog="l.gross@uq.edu.au, version Jan 2021")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration (no py extension)')
#parser.add_argument('--optimize', '-o',  dest='optimize', action='store_true', default=False, help="Calibrated the value for sigma0 before iteration starts.(ignored)")
#parser.add_argument('--vtk', '-v',  dest='vtk', action='store_true', default=False, help="VTK format is used for output otherwise silo is used.")
#parser.add_argument('--outfile', '-o', dest='OUTFILE', metavar='OUTFILE', type=str, default='result', help='name of silo/vtk output file')
#parser.add_argument('--w2', '-2', dest='W2', metavar='W2', type=float, default=1., help='H2 weighting')
parser.add_argument('--gradT', '-g', dest='GradTemp', metavar='GradTemp', type=float, default=30., help='temperature gradient [K/km] at bottom of domain')
parser.add_argument('--tempSurf', '-s', dest='TSURF', metavar='TSURF', type=float, default=15, help='surface temperature, used if no surface temperature file is given')
#parser.add_argument('--extrapolateTemp', '-E',  dest='ExtrapolateTemp', action='store_true', default=False, help="use extrapolated temperatures")
#parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=None, default="restart", help='reststart file name')
#parser.add_argument('--tags_top', '-T', dest='tags_top', metavar='tags_top', type=str, default="top", help='tag of surface elements')
#parser.add_argument('--restart', '-r',  dest='restart', action='store_true', default=False, help="start from restart file. RESTARTFN need to be set and exist.")
parser.add_argument('--debug', '-d',  dest='debug', action='store_true', default=False, help="shows more information.")
args = parser.parse_args()
config = importlib.import_module(args.config)

logger=logging.getLogger('inv')
if args.debug:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

if getMPIRankWorld() == 0:
    print("** This is an Temperature inversion @ %s **"%datetime.now().strftime("%d.%m.%Y %H:%M"))
    print("configuration "+args.config+" imported.")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
if getMPIRankWorld() == 0:
    print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

domain=ReadMesh(config.meshfile)
if getMPIRankWorld() == 0:
    print("Mesh read from "+config.meshfile)


survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns,
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter, commend='#', printInfo=args.debug)
assert survey.getNumObservations()>0, "no data found."

#if getMPIRankWorld() == 0:
#    print("INPUT: w2, output = %e, %s"%(args.W2, args.OUTFILE))

#tags_top=[s.strip() for s in args.tags_top.split(',')  ]
#usetags_top=tags_top is [s.strip()  for s in domain.showTagNames().split(',') if s.strip() in tags_top ]

z = domain.getX()[2]
# create a mask of nodes on the surface:
if config.topsurfaces :
    maskSurfaceTemperatureNodes=MaskFromBoundaryTag(domain, *tuple(config.topsurfaces))
    maskTopSurfaceElements=Scalar(0., FunctionOnBoundary(domain))
    [ maskTopSurfaceElements.setTaggedValue(t, 1.) for t in config.topsurfaces ]
    if getMPIRankWorld() == 0:
        print("surface tag "+str(config.topsurfaces)+" are used.")
else:
    maskSurfaceTemperatureNodes=whereZero(z-sup(z))
    maskTopSurfaceElements=whereZero(FunctionOnBoundary(domain).getX()[2]-sup(z))
# if a surface temperature  file is given the  T_surf, maskTemperatureInterpolatedFromData is extracted from the    
if config.surfacetemperature_file:
    DEMundef=config.surfacetemperature_undef
    dem=np.loadtxt(config.surfacetemperature_file, comments='#', delimiter=',' , skiprows = config.surfacetemperature_skiprows )
    DEMidx=np.logical_not(dem[:,4]==DEMundef)
    DEMx=dem[DEMidx ,1]
    DEMy=dem[DEMidx ,2]
    surfTemperatureData=dem[DEMidx ,4]
    x=domain.getX()
    interp=InterpolatorWithExtension((DEMx, DEMy), surfTemperatureData  , inf(x[0]), sup(x[0]), inf(x[1]), sup(x[1]),  markExtrapolation=False)
    T_surf, maskTemperatureInterpolatedFromData = mapToDomain2D(Scalar(0., maskSurfaceTemperatureNodes.getFunctionSpace()), interp, where=maskSurfaceTemperatureNodes)
    print("surface temperature loaded from ", config.surfacetemperature_file)
else:
    T_surf=Scalar(args.TSURF, Solution(domain))
    maskTemperatureInterpolatedFromData = maskSurfaceTemperatureNodes
    print("surface temperature set to ",str(args.TSURF))

#  temp_set = nodes where temperature is known from surface temperaturs and are ot the top surface
temp_set = maskSurfaceTemperatureNodes * maskTemperatureInterpolatedFromData

T_bg=createBackgroundTemperature(domain, T=T_surf, mask_fixed_temperature=temp_set,
                                 gradT=config.surfacetemperature_gradient,
                                 mask_top_surface=maskTopSurfaceElements)

print(" background temperature = ", str(T_bg))
# ..................
survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=False, columns=[ 'R', 'ETA'],
                     dipoleInjections=True, dipoleMeasurements=True, delimiter=',', commend='#', printInfo=args.debug, unDefined=-999999)
assert survey.getNumObservations()>0, "no data found."
print("data read from "+config.datafile)
print(survey.getNumObservations()," records found.")

# ....create access to the conductivity model ... :
my_conductivity_model=ConductivityModelByTemperature()
print(f"Conductivity model initiated:")
print(f"\tT_ref = {my_conductivity_model.T_ref}.")
print(f"\tsigma_0_ref = {my_conductivity_model.sigma_0_ref}")
print(f"\tP_0 = {my_conductivity_model.P_0}")
print(f"\tsigma_0(1) = {my_conductivity_model.sigma_0_ref/my_conductivity_model.T_ref**my_conductivity_model.P_0}")
print(f"\tMn_ref = {my_conductivity_model.Mn_ref}")
print(f"\tP_Mn = {my_conductivity_model.P_Mn}")
print(f"\tMn(1) = {my_conductivity_model.Mn_ref/my_conductivity_model.T_ref**my_conductivity_model.P_Mn}")

#maskTopSurfaceElements.expand()
#saveSilo("x", n=maskSurfaceTemperatureNodes, e=maskTopSurfaceElements, Tsurf=T_surf, mm=maskTemperatureInterpolatedFromData, Tbg=T_bg)

costf=TemperatureBasedIPInversion(domain=domain, data=survey,


                                  weighting_misfit_0=0.5, pde_tol=1.e-10, stationsFMT="s%s",
                                  w2=args.W2, conductivity_model=my_conductivity_model,
                                  background_temperature=T_bg, mask_fixed_temperature=temp_set, logclip=6, logger=logger)
costf=IPInversion(domain, data=survey,
                  sigma_0_ref=sigma_ref, Mn_ref=Mn_ref, sigma_background=config.sigma_background,
                  w1=config.w1, theta=config.theta, usel1=config.usel1, epsl1=config.epsl1,
                  mask_fixed_property=fixedm,
                  weighting_misfit_ERT=config.weighting_misfit_ERT, pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                  logclip=config.clip_property_function,
                  EPSILON=1e-8, logger=logger)
#saveSilo("t", T_surf=T_surf, temp_set=maskTemperatureInterpolatedFromData, T_bg=T_bg, gradT=grad(T_bg), surface_temperature_weighting=surface_temperature_weighting)
1/0

# this is just for testing:
if False:
        
        x=domain.getX()
        #f=(sup(x[2])-x[2])**2*(x[2]-inf(x[2]))**2
        f=(x[2]-inf(x[2]))**2*(sup(x[1])-x[1])**2*(x[0]-inf(x[0]))**2*(inf(x[2])-x[2])**2
        f/=sup(f)
        c=[ (sup(x[0])+inf(x[0]))/2, (sup(x[1])+inf(x[1]))/2, (sup(x[2])+inf(x[2]))/2]
        L=Lsup(length(domain.getX()))

        m=f
        #m=Scalar(0.2, Solution(domain) ) 
        dd=f*(length(x-c)/(L*0.1))**2
        print("dd,m =",dd,m)

        # test misfit
        #T_surf=15
        #T_min=10.6
        #T_bottom=100.
        #z=domain.getX()[2]
        #a=log((T_bottom-T_min)/(T_surf-T_min))
        #m=a*(z-inf(z))/(sup(z)-inf(z))
        #im=interpolate(m, Function(domain))
        #T=(T_surf-T_min)*exp(im)+T_min
        #print("***** T temperature is reset!")     


        dm=dd
        args=costf.getArgumentsAndCount(m)
        J0=costf.getValueAndCount(m,*args)
        print("J0=",J0)
        
        G=costf.getGradientAndCount(m,*args)
        Dex=costf.getDualProduct(dm,G)

        #print("rel. incr",Lsup(dm)/Lsup(m))
        for k in range(3, 19):
            a=0.5**k
            J1=costf.getValueAndCount(m+a*dm)
            D=(J1-J0)/a
            print("XX", a, "->", J1,J0,  D, Dex, (D-Dex), '->', (D-Dex)/a,)
        1/0


def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, norm_gradFm, args_m, failed):
    if args.RESTARTFN and iterCount >0:
        m.dump(args.RESTARTFN)
        if getMPIRankWorld() == 0: 
            print("restart file %s for step %s created."%(iterCount, args.RESTARTFN))
    #print(f"snapshot for step {k} saved.")
    #saveSilo("snapshot_"+args.OUTFILE, m=x)


# set up solver:

solver= MinimizerLBFGS(F=costf, iterMax=IMAX, logger=logger)
solver.getLineSearch().setOptions(interpolationOrder=INTERPOLATION_ORDER)
solver.setOptions(m_tol=TOL, truncation=TRUNCATION, restart=RESTART, grad_tol=GTOL)
print(solver.getOptions())
#solver.setOptions(interpolationOrder=INTERPOLATION_ORDER, , max_zoom_steps=20,  )#initialHessian=1e9*args.W1)
solver.setCallback(myCallback)
# run solver:
m_init=Scalar(0., Solution(domain) )

# j_tol=TOL
# if args.restart and args.restart:
#   kk=[v for i,v in enumerate(os.listdir()) if v.startswith(args.RESTARTFN) ]
#   if kk:
#     m_init=load(args.RESTARTFN, domain)
#     txt=str(m_init)
#     if getMPIRankWorld() == 0:
#             print("restart file %s read. initial M = %s"%(args.RESTARTFN, txt))
#     j_tol=None
#
# solver.setTolerance(m_tol=TOL, J_tol=j_tol)

solver.run(m_init)
solver.logSummary()
m=solver.getResult()

T=costf.getTemperature(m)
    
txt4=str(T)
if getMPIRankWorld() == 0: 
        print("final T =", txt4)

 
 
outargs={ "tag" : makeTagField(Function(domain)), "T" : T, "S" : m }
if args.vtk:
    saveVTK(args.OUTFILE, **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.vtu"%args.OUTFILE)
else:
    saveSilo(args.OUTFILE,  **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.silo"%args.OUTFILE)

if getMPIRankWorld() == 0: 
    print("All done - Have a nice day!")
