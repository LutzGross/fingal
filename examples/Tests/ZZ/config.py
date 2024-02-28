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
datacolumns = ['R', 'ETA']
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
eta_ref= 0.01
Mn_ref=eta_ref/(1-eta_ref)*sigma_ref
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
m_tolerance=1.e-3
g_tolerance=1.e-6
interpolation_order=3
imax=400
truncation=20
restart=60
pde_tol=1e-10
w1=1E-3
useL1Norm=False
epsilonL1Norm=1e-4
useLogMisfit = False

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


