#
# this is a fingal configuration file
#
created = 'today'
project='Wanner'
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
sigma_ref=0.02
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
g_tolerance=1.e-4
interpolation_order=3
imax=400
truncation=20
restart=60
pde_tol=1e-10

useL1Norm=False
epsilonL1Norm=1e-4

w1=1.e-4
w1=1e-6 # ok.
w1=1e-5 # not bad
w1=1e-5 # ok.


useLogMisfit = True

#w1=5.e-7
#useLogMisfit = False


#
# Output handeling:
#
outfile='sigma'
core = ['core']
#   tag(s) for face elements (excluding top surface)
faces = ['faces']
restartfile = 'restart'


