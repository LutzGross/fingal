#!/usr/bin/python3
import sys

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.Get_size() < 1 :
        raise ValueError("This runs on single node only.")
except ImportError:
    pass
from esys.escript import *
from esys.finley import ReadGmsh, ReadMesh
from fingal import readElectrodeLocations, readSurveyData, ConductivityModelByTemperature, IPSynthetic, makeMaskForOuterSurface, setupERTPDE
from esys.weipa import saveVTK, saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import getMaskFromBoundaryTag, Locator

import numpy as np
import argparse
import importlib
import subprocess


parser = argparse.ArgumentParser(description='generates a data file from temperature distribution.', epilog="version 3/2024")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--coredepth', '-d',  dest='coredepth', type=int, default=35, help="core depth relative to core width around core in %% of electrode area size (default 35)")
parser.add_argument('--extracore', '-e',  dest='extracore', type=int, default=20, help="relative extracore padding within the core around electrodes in %%  of edge legth (default 20)")
parser.add_argument('--padding', '-p',  dest='padding', type=int, default=1100, help="relative padding around core in %% of core edge length (default 1100)")
parser.add_argument('--coremeshfactor', '-C',  dest='coremeshfactor', type=float, default=0.5, help="refinement factor of mesh in core relative to the electrode distance (default 0.5)")
parser.add_argument('--stationmeshfactor', '-s',  dest='stationmeshfactor', type=float, default=0.3, help="refinement factor of mesh near stations relative to their distance (default 0.3)")
parser.add_argument('--paddingmesh', '-P',  dest='paddingmesh', type=float, default=20, help="number of element on the longest edge of the padding region (default 20)")
parser.add_argument('--geofile', '-g',  dest='geofile', type=str, default="truedomain", help="name of gmsh geofile to generate")
parser.add_argument('--haveMesh', '-M',  dest='haveMesh', action='store_true', default=False, help="if set mesh generation is not started. Useful for debugging.")
parser.add_argument('--silofile', '-o', dest='silofile', metavar='SILOFILE', type=str, default='truemodel', help='name of silo output file (generated)')
parser.add_argument('--noise', '-n', dest='noise', metavar='NOISE', type=float, default=0, help='%% of noise added to data.')
#parser.add_argument('--conductiveSurface', '-c',  dest='conductiveSurface', action='store_true', default=False, help="use conductive surface. Otherwise surface temperature is fixed.")

args = parser.parse_args()
config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

MSHFILE=args.geofile + ".msh"
FLYFILE=args.geofile + ".fly"
GEOFILE = args.geofile + ".geo"
#=====================
T_air=15  # air temperature in celcius!
q_bottom= config.background_heat_surface_flux_bottom
h_top=0.05*2 # convective heat transfer coefficient in  W/m/K
Q_volcano = 0.8# heatsource of volcane W/m^3/K
#===========================================


elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

positions=np.array([ [ elocations[s][0], elocations[s][1], elocations[s][2]] for s in elocations ] )
zCore=np.mean(positions[:,2])
zmin=positions[:,2].min()
xmin=positions[:,0].min()
xmax=positions[:,0].max()
ymin=positions[:,1].min()
ymax=positions[:,1].max()


DistanceElectrodes=1e999
for i in range(len(elocations)):
    for j in range(len(elocations)):
        if i < j:
            d=np.linalg.norm(positions[j]-positions[i])
            DistanceElectrodes=min(DistanceElectrodes, d)

print("level of survey = ", zCore)
print("x range = ", xmin, xmax)
print("y range = ", ymin, ymax)
print("DistanceElectrodes = ", DistanceElectrodes)

if not args.haveMesh:

    diagonalAreaOfElectrodes = ((xmax - xmin) ** 2 + (ymax - ymin) ** 2) ** 0.5
    assert diagonalAreaOfElectrodes > 0., "area of electrodes is zero."
    XextensionCore = (max(xmax - xmin, ymax - ymin) * args.extracore) / 100.
    YextensionCore = XextensionCore

    if xmax - xmin < 1e-5 * diagonalAreaOfElectrodes:
        XextensionCore = YextensionCore
    if ymax - ymin < 1e-5 * diagonalAreaOfElectrodes:
        YextensionCore = XextensionCore

    CoreThickness = (diagonalAreaOfElectrodes * args.coredepth) / 100.
    Xpadding = ((xmax - xmin + 2 * XextensionCore) * args.padding) / 100.
    Ypadding = ((ymax - ymin + 2 * YextensionCore) * args.padding) / 100.
    Zpadding = min(Xpadding, Ypadding)
    Xpadding, Ypadding = Zpadding, Zpadding
    XminBB = xmin - XextensionCore - Xpadding
    XmaxBB = xmax + XextensionCore + Xpadding
    YminBB = ymin - YextensionCore - Ypadding
    YmaxBB = ymax + YextensionCore + Ypadding
    ZminBB = zCore - CoreThickness - Zpadding
    ZminCore = zCore - CoreThickness

    out = "// Core:\n"
    out += "XminCore = %g;\n" % (xmin - XextensionCore)
    out += "XmaxCore = %g;\n" % (xmax + XextensionCore)
    out += "YminCore = %g;\n" % (ymin - YextensionCore)
    out += "YmaxCore = %g;\n" % (ymax + YextensionCore)
    out += "ZCore = %g;\n" % zCore
    out += "ZminCore = %g;\n" % ZminCore
    out += "CoreThickness = %g;\n" % (zCore-ZminCore)
    out += "// Big Box\n"
    out += "XminBB = %g;\n" % XminBB
    out += "XmaxBB = %g;\n" % XmaxBB
    out += "YminBB = %g;\n" % YminBB
    out += "YmaxBB = %g;\n" % YmaxBB
    out += "DepthBoundingBox = %g;\n" % (CoreThickness + Zpadding)
    out += "// Geometry:\n"
    out += "DistanceElectrodes = %g;\n" % (DistanceElectrodes)
    out += "DepthInterface = 3 * DistanceElectrodes;\n"
    out += "RadiusConduit = 3 * DistanceElectrodes;\n"
    out += "OffsetXConduit = 3 * DistanceElectrodes;\n"
    out += "OffsetYConduit = 0.0 * DistanceElectrodes;\n"
    out += "// mesh sizes\n"
    out += "meshSizeCore = %g * DistanceElectrodes;\n" % (args.coremeshfactor)
    out += "meshSizeBB = %g/%d;\n" % (max(xmax - xmin + 2 * XextensionCore + 2 * Xpadding, ymax - ymin + 2 * YextensionCore + 2 * Ypadding),
    args.paddingmesh)
    out += "meshSizeElectrodes = DistanceElectrodes * %s;\n" % (args.stationmeshfactor)
    out += "meshSizeConduit = meshSizeCore;\n"

    GEOTEMPLATE=open("truegeo_template.geo", 'r').read()

    ELECTRODES = ""

    for j, e  in enumerate(elocations):
        ID=config.stationsFMT%e
        ELECTRODES+= f"Point(k + {j+1}) = {{{elocations[e][0]}, {elocations[e][1]}, {elocations[e][2]}, meshSizeElectrodes}};\n"
        ELECTRODES+= f"Point{{k + {j+1}}} In Surface {{355}};\n"
        ELECTRODES+= f'Physical Point("{ID}")  = {{ k+{j+1} }};\n'

    open(GEOFILE, 'w').write(GEOTEMPLATE.replace('%ELECTRODES%', ELECTRODES).replace('%HEADER%', out))
    print("true geometry file written to "+GEOFILE)


    rp=subprocess.run(["gmsh", "-3", "-optimize_netgen", "-o", MSHFILE, GEOFILE])
    #rp.check_returncode()
    print(">> GMSH mesh file %s was generated."%MSHFILE)



    # read mesh
    dts=[]
    dps=[]
    for s in elocations:
        if config.stationsFMT:
            dts.append(config.stationsFMT%s)
        else:
            dts.append(s)
        dps.append(elocations[s])
    domain=ReadGmsh(MSHFILE, 3, diracPoints=dps, diracTags=dts, optimize=True )
    domain.write(FLYFILE)
else:
    domain=ReadMesh(FLYFILE)
    print("mesh read from "+FLYFILE)
#---------- construct temperature --------------------------------
pde = setupERTPDE(domain, tolerance=1e-10, poisson=False, debug=1)
pde.setValue(A=config.thermal_conductivity *kronecker(3))

# flux at top and buttomqs
h=Scalar(0., FunctionOnBoundary(domain))
h.setTaggedValue("SurfaceTop", h_top)
qs= h * T_air
qs.setTaggedValue("SurfaceBottom", -q_bottom)
pde.setValue(y=qs, d=h)

Q=Scalar(config.background_heat_production, Function(domain))
Q.setTaggedValue("Volcano", Q_volcano)
pde.setValue(Y=Q)

T=pde.getSolution()

saveSilo("temperature", T=T, flux =  -config.thermal_conductivity *grad(T))
del pde
# save surface temperature to file:
if False:
    T_surf= 15
    T_bottom = 20
    z=domain.getX()[2]
    #iz=interpolate(z, Function(domain))
    T=(T_surf-T_bottom)/(sup(z)-inf(z))*(z-inf(z))+T_bottom
    print("***** T temperature is reset!") 
#    print(((sup(x)+inf(x))/2)**2)
#    T=(T_bottom-T_air)/(inf(z)-sup(z))*(z-sup(z))
#    T*=(sup(x)-x)*(x-inf(x))/((sup(x)-inf(x))/2)**2
#    T*=(sup(y)-y)*(y-inf(y))/((sup(y)-inf(y))/2)**2AC2069
#    T+=T_air
print("Temperature =",str(T))   

if config.surfacetemperature_file:
    x=domain.getX()

    saveDataCSV(config.surfacetemperature_file, A0=x[0], A1=x[1], A3=x[2], T=T, mask=whereZero(x[2]-sup(x[2])), refid=True)
    print("surface temperature saved to "+config.surfacetemperature_file)

if config.surfaceflux_file:
    qq = h * ( T - T_air)
    mm = Scalar(0., qq.getFunctionSpace())
    mm.setTaggedValue("SurfaceTop", 1)
    x=qq.getX()

    saveDataCSV(config.surfaceflux_file, A0=x[0], A1=x[1], A3=x[2], Qr=qq, mask=mm, refid=True)
    print("surface flux saved to "+config.surfaceflux_file)


# get the schedule
survey=readSurveyData(config.schedulefile, stations=elocations,
                      usesStationCoordinates=config.usesStationCoordinates, columns=[],
                      dipoleInjections=config.dipoleInjections,
                      dipoleMeasurements=config.dipoleMeasurements,
                      delimiter=config.datadelimiter,
                      commend='#', printInfo=True)

# set conductivity:
conductivity_model=ConductivityModelByTemperature()

iT=interpolate(T, Function(domain))
sigma0=conductivity_model.getDCConductivity(T=iT)
Mn=conductivity_model.getChargeability(T=iT)

loc=Locator(T.getFunctionSpace(),  [ survey.getStationLocation(s) for s in survey.getStationNumeration()] )
iT_at_stations = np.array(loc(T))
sigma0_at_stations=conductivity_model.getDCConductivity(T=iT_at_stations)
Mn_at_stations=conductivity_model.getChargeability(T=iT_at_stations)


print("sigma0 = ",str(sigma0))
print("sigma0_at_stations = ",str(sigma0_at_stations))
print("Mn = ",str(Mn))
print("Mn_at_stations = ",str(Mn_at_stations))
# -----------------------------------------------------------------------------------
tag_faces = ['SurfaceBottom', 'VerticalFaces']

mask_zero_potential = getMaskFromBoundaryTag(domain, *tuple(tag_faces) )

runner=IPSynthetic(domain, survey,  sigma_src=conductivity_model.sigma_0_ref,
                    maskZeroPotential = mask_zero_potential,
                    stationsFMT=config.stationsFMT,
                    createSecondaryData=True,
                    createFieldData=False,  printInfo = True)

runner.setProperties(sigma_0=sigma0, M_n=Mn,
                     sigma_0_at_stations=sigma0_at_stations, M_n_at_stations = Mn_at_stations)

runner.write(config.datafile , datacolumns = config.datacolumns, addNoise = args.noise>0,
                            rel_error = args.noise, delimiter=config.datadelimiter,
                            usesStationCoordinates= config.usesStationCoordinates,
                            iFMT="%d", dFMT="%.5g", xFMT="%e")
# -----------------------------------------------------------------------------------
if args.silofile:
    kwargs = { "Mn" : Mn, "sigma0" : sigma0, "eta" : Mn/(Mn+sigma0), "tag" : getRegionTags(Function(domain)) }
    A=list(elocations.keys())[0]
    iA=survey.getStationNumber(A)
    kwargs[f'V0{A}'] = runner.potential_0[iA]
    kwargs[f'V2{A}'] = runner.potential_0[iA]-runner.potential_oo[iA]
    kwargs[f'ETA{A}'] = kwargs[f'V2{A}']/kwargs[f'V0{A}']
    if runner.createFieldData:
        kwargs[f'E0_s{A}'] = runner.field_0[iA] + runner.source_field[iA]
        kwargs[f'E2_s{A}'] = runner.field_0[iA] - runner.field_oo[iA]
    saveSilo(args.silofile , **kwargs)
    print(f'values {kwargs.keys()} written to file {args.silofile}.')

