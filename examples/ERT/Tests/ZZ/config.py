#
# this is a fingal configuration file
#
created = '31.08.2023 15:10'
project='ZZ'
# 
#   name of the mesh file
# 
meshfile = 'mesh.fly'
#
#  file of the location of stations/electrodes. 
#
stationfile = 'stations.csv'
stationdelimiter = ','
stationsFMT = 's%s'
L_stations=1.
adjustStationLocationsToElementCenter=True
#
# this defines the data file: 
# 
datafile = 'data.csv'
datacolumns = ['R']
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = 'schedule.csv'
#
#  This section of the file defines the inversion
#

region_fixed=[]
sigma_ref=0.005
Mn_ref=0.01/(1-0.01)*sigma_ref
#def true_properties(domain):
#    from esys.escript import Scalar, Function
#    sigma_true=Scalar(sigma_ref, Function(domain))
#    Mn_true=Scalar(Mn_ref, Function(domain))
#    return sigma_true, Mn_true
#
#
#  Inversion:
#
weighting_misfit_ERT=0.5
clip_property_function=10
m_tolerance=1.e-4
g_tolerance=1.e-4
interpolation_order=3
imax=400
truncation=20
restart=60
pde_tol=1e-10
w1=1e-4
usel1=False
epsl1=1e-4
w0=0.0
theta=0
alpha0=1
alpha1=1.0
w1gamma=w1
w0gamma=w0
alpha0gamma=alpha0
alpha1gamma=alpha1
#
# Output handeling:
#
outfile='sigma'
core = ['core']
#   tag(s) for face elements (excluding top surface)
faces = ['faces']
restartfile = 'restart'

def true_properties(domain):
    DD=5 # depth of high conductivity area.
    KAPPA=3 # curvture factor
    x0=30 # horizontal offset of high conductivit area
    InterfaceDepth=35
    sigma_low =sigma_ref
    sigma_high = sigma_ref*100
    from esys.escript import Scalar, ReducedFunction, interpolate, Function, whereNegative, wherePositive
    X=ReducedFunction(domain).getX()
    z=X[2]
    x=X[0]
    h=-DD-KAPPA*(x-x0)**2
    mp=wherePositive(whereNegative(z-h)+whereNegative(z+InterfaceDepth))
    sigma_true=sigma_high * mp + (1-mp) * sigma_low
    Mn_true=Scalar(Mn_ref, Function(domain))
    return interpolate(sigma_true, Function(domain)), Mn_true
