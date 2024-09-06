#!/usr/bin/python3
from esys.escript import *
from esys.finley import ReadGmsh
from fingal import readElectrodeLocations, readSurveyData, ConductivityModelByTemperature, IPSynthetic, makeMaskForOuterSurface
from esys.weipa import saveVTK, saveSilo
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import MaskFromTag, Locator
import numpy as np
import argparse
import importlib
import subprocess


parser = argparse.ArgumentParser(description='generates a data file from temperature distribution.', epilog="version 3/2024")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--coredepth', '-d',  dest='coredepth', type=int, default=45, help="core depth relative to core width around core in %% of electrode area size (default 45)")
parser.add_argument('--extracore', '-e',  dest='extracore', type=int, default=20, help="relative extracore padding within the core around electrodes in %%  of edge legth (default 20)")
parser.add_argument('--padding', '-p',  dest='padding', type=int, default=150, help="relative padding around core in %% of core eddge length (default 150)")
parser.add_argument('--coremeshfactor', '-C',  dest='coremeshfactor', type=float, default=1., help="refinement factor of mesh in core relative to the electrode distance (default 0.1)")
parser.add_argument('--stationmeshfactor', '-s',  dest='stationmeshfactor', type=float, default=0.3, help="refinement factor of mesh near stations relative to their distance (default 0.3)")
parser.add_argument('--paddingmesh', '-P',  dest='paddingmesh', type=float, default=20, help="number of element on the longest edge of the padding region (default 20)")
parser.add_argument('--geofile', '-g',  dest='geofile', type=str, default="truedomain", help="name of gmsh geofile to generate")
parser.add_argument('--noTrueMesh', '-N',  dest='noTrueMesh', action='store_true', default=False, help="if set only gmsh geo file is generated but mesh generation is not started. Useful for debugging.")
parser.add_argument('--silofile', '-o', dest='silofile', metavar='SILOFILE', type=str, default='truemodel', help='name of silo output file (generated)')
parser.add_argument('--noise', '-n', dest='noise', metavar='NOISE', type=float, default=0, help='%% of noise added to data.')
#parser.add_argument('--conductiveSurface', '-c',  dest='conductiveSurface', action='store_true', default=False, help="use conductive surface. Otherwise surface temperature is fixed.")

args = parser.parse_args()
config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")

#=====================
T_air=15  # air temperature in celcius!
T_volcano=150 # temperature in volcano\
K_top, K_base=1.,1.
Q_bottom=config.surfacetemperature_gradient * K_top
h_top=0.1 # convective heat transfer coefficient in  W/mK
#===========================================

MSHFILE=args.geofile+".msh"


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

if not args.noTrueMesh:

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
    out += "meshSizeConduit = meshSizeCore /2 ;\n"

    GEOTEMPLATE=open("truegeo_template.geo", 'r').read()

    ELECTRODES = ""

    for j, e  in enumerate(elocations):
        ID=config.stationsFMT%e
        ELECTRODES+= f"Point(k + {j+1}) = {{{elocations[e][0]}, {elocations[e][1]}, {elocations[e][2]}, meshSizeElectrodes}};\n"
        ELECTRODES+= f"Point{{k + {j+1}}} In Surface {{355}};\n"
        ELECTRODES+= f'Physical Point("{ID}")  = {{ k+{j+1} }};\n'
    GEOFILE=args.geofile+".geo"
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
#---------- construct temperature --------------------------------
pde=LinearSinglePDE(domain, isComplex=False)
pde.setSymmetryOn()
optionsG=pde.getSolverOptions()   
optionsG.setTolerance(1e-10)



# thermal conductivity:
K=Scalar(K_base, Function(domain))
K.setTaggedValue("LayerTop", K_top)
K.setTaggedValue("LayerBottom", K_base)
pde.setValue(A=K*kronecker(3))
# flux at top and buttom
Q=Scalar(0., FunctionOnBoundary(domain))
Q.setTaggedValue("SurfaceBottom", Q_bottom)
Q.setTaggedValue("SurfaceTop", h_top*T_air)
h=Scalar(0., FunctionOnBoundary(domain))
h.setTaggedValue("SurfaceTop", h_top)
Q.expand()
h.expand()
pde.setValue(y=Q, d=h)
q=MaskFromTag(domain, "Volcano")
pde.setValue(q=q, r=T_volcano*q)
T=pde.getSolution()
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
#    T*=(sup(y)-y)*(y-inf(y))/((sup(y)-inf(y))/2)**2
#    T+=T_air
print("Temperature =",str(T))   

if config.surfacetemperature_file:
    x=domain.getX()

    saveDataCSV(config.surfacetemperature_file, A0=x[0], A1=x[1], A3=x[2], T=T, mask=whereZero(x[2]-sup(x[2])), refid=True)
    print("surface temperature saved to "+config.surfacetemperature_file)


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

iT_faces=interpolate(T, FunctionOnBoundary(domain))
sigma0_faces=conductivity_model.getDCConductivity(T=iT_faces)
Mn_faces=conductivity_model.getChargeability(T=iT_faces)

loc=Locator(T.getFunctionSpace(),  [ survey.getStationLocation(s) for s in survey.getStationNumeration()] )
iT_at_stations = np.array(loc(T))
sigma0_at_stations=conductivity_model.getDCConductivity(T=iT_at_stations)
Mn_at_stations=conductivity_model.getChargeability(T=iT_at_stations)


print("sigma0 = ",str(sigma0))
print("sigma0_faces = ",str(sigma0_faces))
print("Mn_faces = ",str(Mn_faces))
print("Mn = ",str(Mn))

# -----------------------------------------------------------------------------------
tag_faces = ['SurfaceBottom', 'VerticalFaces']
runner=IPSynthetic(domain, survey,  sigma_src=config.sigma_ref,
                    mask_faces = makeMaskForOuterSurface(domain, taglist=tag_faces),
                    stationsFMT=config.stationsFMT,
                    createSecondaryData=True,
                    createFieldData=False,  printInfo = True)
runner.mask_faces.expand()

runner.setProperties(sigma_0=sigma0, sigma_0_faces=sigma0_faces,
                     M_n=Mn, M_n_faces=Mn_faces,
                     sigma_0_at_stations=sigma0_at_stations, M_n_at_stations = Mn_at_stations)

runner.write(config.datafile , datacolumns = config.datacolumns, addNoise = args.noise>0,
                            rel_error = args.noise, delimiter=config.datadelimiter,
                            usesStationCoordinates= config.usesStationCoordinates,
                            iFMT="%d", dFMT="%.5g", xFMT="%e")
# -----------------------------------------------------------------------------------
if args.silofile:
    kwargs = { "Mn" : Mn, "sigma0" : sigma0, "tag" : makeTagMap(Function(domain)) }
    A=list(elocations.keys())[0]
    iA=survey.getStationNumber(A)
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

