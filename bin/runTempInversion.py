#!/usr/bin/python3
try:
    from mpi4py import MPI
except ImportError:
    pass
import importlib, sys, os

from esys.escript import *


from datetime import datetime
from esys.escript.minimizer import MinimizerLBFGS
sys.path.append(os.getcwd())
from fingal import readElectrodeLocations, readSurveyData, InterpolatorWithExtension, mapToDomain2D, createBackgroundTemperature
from fingal import ConductivityModelByTemperature, InversionIPBySource, InversionIPByFlux, InversionIPByTemperature, InversionIPByFluxWithWeakBC
from esys.finley import ReadMesh
import numpy as np
from esys.weipa import saveSilo
import argparse
from esys.escript.pdetools import getMaskFromBoundaryTag
import logging
#from IPtemperatureinversion import ConductivityModelByTemperature, createBackgroundTemperature, TemperatureBasedIPInversion, InterpolatorWithExtension, mapToDomain2D


parser = argparse.ArgumentParser(description='driver to invert an electric field intensity survey (aka Fullwaver). Measurements may not be dipoles.', epilog="l.gross@uq.edu.au, version Jan 2021")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration (no py extension)')
#parser.add_argument('--optimize', '-o',  dest='optimize', action='store_true', default=False, help="Calibrated the value for sigma0 before iteration starts.(ignored)")
parser.add_argument('--vtk', '-v',  dest='vtk', action='store_true', default=False, help="VTK format is used for output otherwise silo is used.")
#parser.add_argument('--outfile', '-o', dest='OUTFILE', metavar='OUTFILE', type=str, default='result', help='name of silo/vtk output file')
#parser.add_argument('--w2', '-2', dest='W2', metavar='W2', type=float, default=1., help='H2 weighting')
#parser.add_argument('--gradT', '-g', dest='GradTemp', metavar='GradTemp', type=float, default=30., help='temperature gradient [K/km] at bottom of domain')
#parser.add_argument('--extrapolateTemp', '-E',  dest='ExtrapolateTemp', action='store_true', default=False, help="use extrapolated temperatures")
parser.add_argument('--restartfile', '-R', dest='RESTARTFN', metavar='RESTARTFN', type=None, default="restart", help='reststart file name')
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
if getMPIRankWorld() == 0:
    print("data read from "+config.datafile)
    print(survey.getNumObservations()," records found.")


# create a mask of nodes on the surface:
z = domain.getX()[2]
if config.surface_tags :
    maskSurfaceTemperatureNodes=getMaskFromBoundaryTag(domain, *tuple(config.surface_tags))
    maskTopSurfaceElements=Scalar(0., FunctionOnBoundary(domain))
    [ maskTopSurfaceElements.setTaggedValue(t, 1) for t in config.surface_tags]
    if getMPIRankWorld() == 0:
        print("surface tags =  "+str(config.surface_tags)+".")
else:
    maskSurfaceTemperatureNodes=whereZero(z-sup(z))
    maskTopSurfaceElements=wherePositive(interpolate(whereZero(ReducedFunctionOnBoundary(domain).getX()[2] - sup(z)), FunctionOnBoundary(domain)))
maskBottomSurfaceElements=wherePositive(interpolate(whereZero(ReducedFunctionOnBoundary(domain).getX()[2] - inf(z)), FunctionOnBoundary(domain)))
# ====================================================================
# ....create access to the conductivity model ... :
my_conductivity_model=ConductivityModelByTemperature()
if getMPIRankWorld() == 0:
    print(f"Conductivity model initiated: T_ref = {my_conductivity_model.T_ref}\n\t\t\tsigma_0_ref = {my_conductivity_model.sigma_0_ref}, exponent = {my_conductivity_model.P_0}\n\t\t\tMn_ref = {my_conductivity_model.Mn_ref}, exponent = {my_conductivity_model.P_Mn}")


# if a surface temperature  file is given the  T_surf, maskTemperatureInterpolatedFromData is extracted from the    
if config.surfacetemperature_file:
    TEMPundef=config.surfacetemperature_undef
    TEMP=np.loadtxt(config.surfacetemperature_file, comments='#', delimiter=',' , skiprows = config.surfacetemperature_skiprows )
    TEMPidx=np.logical_not(TEMP[:,4]==TEMPundef)
    TEMPx=TEMP[TEMPidx ,1]
    TEMPy=TEMP[TEMPidx ,2]
    surfTemperatureData=TEMP[TEMPidx ,4]
    x=domain.getX()
    interpolator_Temperature=InterpolatorWithExtension((TEMPx, TEMPy), surfTemperatureData  , inf(x[0]), sup(x[0]), inf(x[1]), sup(x[1]),  markExtrapolation=not config.useTemperatureExtrapolation)
    T_surf, maskTemperatureInterpolatedFromData = mapToDomain2D(Scalar(0., maskSurfaceTemperatureNodes.getFunctionSpace()), interpolator_Temperature, where=maskSurfaceTemperatureNodes)
    if getMPIRankWorld() == 0:
        print("surface temperature loaded from ", config.surfacetemperature_file)
    del TEMP, interpolator_Temperature
else:
    T_surf=Scalar(my_conductivity_model.T_ref, Solution(domain))
    maskTemperatureInterpolatedFromData = maskSurfaceTemperatureNodes
    if getMPIRankWorld() == 0:
        print("surface temperature set to ",str(my_conductivity_model.T_ref))
#
#... set conductivity, flux and heat source for background temperature:
#
thermal_conductivity = Scalar(config.thermal_conductivity, Function(domain))
Q = Scalar(config.background_heat_production, Function(domain))
[Q.setTaggedValue(t, config.background_heat_production_core) for t in config.core_tags]
logger.debug("Heat source set to " + str(Q) )
q = - maskBottomSurfaceElements * config.background_heat_surface_flux_bottom + maskTopSurfaceElements * config.background_heat_surface_flux_top
#q = - maskBottomSurfaceElements * config.background_heat_surface_flux_bottom + maskTopSurfaceElements * config.background_heat_surface_flux_top

## read surface fluxes and merge into q:
if config.surfaceflux_file:
    FLUXundef=config.surfaceflux_undef
    FLUX=np.loadtxt(config.surfaceflux_file, comments='#', delimiter=',' , skiprows = config.surfaceflux_skiprows )
    FLUXidx=np.logical_not(FLUX[:,4]==FLUXundef)
    FLUXx=FLUX[FLUXidx ,1]
    FLUXy=FLUX[FLUXidx ,2]
    surfFLUXData=FLUX[FLUXidx ,4]
    x=FunctionOnBoundary(domain).getX()
    interpolator_Flux=InterpolatorWithExtension((FLUXx, FLUXy), surfFLUXData  , inf(x[0]), sup(x[0]), inf(x[1]), sup(x[1]),  markExtrapolation=not config.useFluxExtrapolation)
    Flux_surf, maskFluxInterpolatedFromData = mapToDomain2D(Scalar(0., x.getFunctionSpace()), interpolator_Flux, where=maskTopSurfaceElements)
    if getMPIRankWorld() == 0:
        print("surface flux loaded from ", config.surfaceflux_file)
    logger.debug("Read surface heat fluxes are " + str(Flux_surf))
    q =  maskFluxInterpolatedFromData * Flux_surf + (1-maskFluxInterpolatedFromData) * q
    del FLUX, interpolator_Flux

logger.debug("Surface flux set to " + str(q) )

# ... create mask where the elect. potential is set to zero:
mask_face = getMaskFromBoundaryTag(domain, *config.faces_tags)

# ------------------------------------------------------------------------------------------------
# ===== Start inversion
if config.model == "SOURCE":
    costf = InversionIPBySource(domain, data=survey,  maskZeroPotential=mask_face,
                                conductivity_model=my_conductivity_model,
                                pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                                useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC = config.data_rtol,
                                useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                                weightingMisfitDC=config.regularization_weighting_DC_misfit,
                                w1=config.regularization_w1DC, reg_tol=None,
                                surface_temperature=T_surf, mask_surface_temperature=maskTemperatureInterpolatedFromData,
                                thermal_conductivity=thermal_conductivity,
                                background_heat_source = Q,
                                surface_flux = q, logclip = config.clip_property_function,
                                logger=logger.getChild(f"IP-Temp-SRC"))
    m_init = Scalar(0., Solution(domain))
elif config.model == "FLUXWEAK":
    costf = InversionIPByFluxWithWeakBC(domain, data=survey,  maskZeroPotential=mask_face,
                                conductivity_model=my_conductivity_model,
                                sigma_src=None,
                                pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, length_scale=config.regularization_length_scale,
                                useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC = config.data_rtol,
                                useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                                weightingMisfitDC=config.regularization_weighting_DC_misfit,
                                w1=config.regularization_w1DC, reg_tol=None,
                                surface_flux=q,
                                set_temperature = config.temperature_longrange, mask_set_temperature = mask_face  * maskSurfaceTemperatureNodes ,
                                mask_flux_variation=maskTopSurfaceElements,
                                flux_variation_factor=config.surface_values_variation_factor,
                                thermal_conductivity=thermal_conductivity,
                                heat_production = Q,
                                logger=logger.getChild(f"IP-Temp-FLXBC"))
    m_init = Vector(0., Solution(domain))

elif config.model == "FLUX":
    costf = InversionIPByFlux(domain, data=survey,  maskZeroPotential=mask_face,
                                conductivity_model=my_conductivity_model,
                                sigma_src=None,
                                pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, length_scale=config.regularization_length_scale,
                                useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC = config.data_rtol,
                                useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                                weightingMisfitDC=config.regularization_weighting_DC_misfit,
                                w1=config.regularization_w1DC, reg_tol=None,
                                surface_temperature=T_surf, mask_surface_temperature=maskTemperatureInterpolatedFromData,
                                thermal_conductivity=thermal_conductivity,
                                heat_production = Q,
                                surface_flux = q,
                                logger=logger.getChild(f"IP-Temp-FLX"))
    m_init = Vector(0., Solution(domain))
elif config.model == "TEMP":
    costf = InversionIPByTemperature(domain, data=survey,  maskZeroPotential=mask_face,
                                conductivity_model=my_conductivity_model, logclip = config.clip_property_function,
                                sigma_src=None,
                                pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                                useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC = config.data_rtol,
                                useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                                weightingMisfitDC=config.regularization_weighting_DC_misfit,
                                w1=config.regularization_w1DC, reg_tol=None, temperature_variation_factor = config.surface_values_variation_factor,
                                surface_temperature=T_surf, mask_surface_temperature=maskTemperatureInterpolatedFromData,
                                thermal_conductivity=thermal_conductivity,
                                heat_production = Q,
                                surface_flux = q,
                                logger=logger.getChild(f"IP-Temp-TEMP"))
    m_init = Scalar(0., Solution(domain))
else:
    raise ValueError("Unknown model %s"%config.model)

def myCallback(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
    saveSilo(f"tmp/out{iterCount:03d}", m = m, T= args_m[0], tag =  getRegionTags(ReducedFunction(m.getDomain()) ))
    if getMPIRankWorld() == 0:
        print(f">>>> iteration saved to  file tmp/out{iterCount:03d}!")
    if args.RESTARTFN and iterCount >0:
        m.dump(args.RESTARTFN)
        if getMPIRankWorld() == 0: 
            print("restart file %s for step %s created."%(iterCount, args.RESTARTFN))

def myCallbackQ(iterCount, m, dm, Fm, grad_Fm, norm_m, args_m, failed):
        saveSilo(f"tmpQ/out{iterCount:03d}", m=m, Q=args_m[0], T=args_m[1], tag=getRegionTags(ReducedFunction(m.getDomain())))
        if getMPIRankWorld() == 0:
            print(f">>>> iteration saved to  file tmpQ/out{iterCount:03d}!")
        # if args.RESTARTFN and iterCount > 0:
        #     m.dump(args.RESTARTFN)
        #     if getMPIRankWorld() == 0:
        #         print("restart file %s for step %s created." % (iterCount, args.RESTARTFN))

    #print(f"snapshot for step {k} saved.")
    #saveSilo("snapshot_"+args.OUTFILE, m=x)


# set up solver:

solver= MinimizerLBFGS(F=costf, iterMax=config.imax, logger=logger)
solver.getLineSearch().setOptions(interpolationOrder=config.interpolation_order, iterMax=50)
solver.setOptions(m_tol=config.m_tolerance, truncation=config.truncation, restart=config.restart, grad_tol=config.g_tolerance)
if config.model == "SOURCE":
    solver.setCallback(myCallbackQ)
else:
    solver.setCallback(myCallback)

# run solver:


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
if config.model == "SOURCE":
    iQ = costf.getHeatSource(m, applyInterpolation=True)
    T = cost.getTemperature(iQ)
    V = - grad(T)
elif config.model == "FLUX":
    V = -costf.extractPropertyFunction(m)
    T = costf.getTemperature(V)
    iQ = -div(V)
else:
    T = costf.getTemperature(m)
    V = - grad(T)
txt4=str(T)
if getMPIRankWorld() == 0: 
        print("final T =", txt4)

outargs={ "tag" : getRegionTags(ReducedFunction(domain)), "T"  : T, "maskSetT" : maskTemperatureInterpolatedFromData ,"T_bg" : costf.T_background,
          "dT_ratio" : T/costf.T_background, "V"  : V }
if config.model in [ "FLUX" , "SOURCE" ] :
    outargs["Q"] = iQ

if args.vtk:
    saveVTK(config.outfile, **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.vtu"%config.outfile)
else:
    saveSilo(config.outfile,  **outargs)
    if getMPIRankWorld() == 0: print("result written to %s.silo"%config.outfile)

if getMPIRankWorld() == 0: 
    print("All done - Have a nice day!")
