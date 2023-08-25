#
# this is a fingal configuration file
#
created = '22.08.2023 16:13'
project='Wenner'
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
schedulefile = 'data.csv'
#
#  This section of the file defines the inversion
#

region_fixed=[]
sigma_background=0.0006*0+0.02132941887325601
sigma_ref=sigma_background
Mn_ref=0.01/(1-0.01)*sigma_background
#def true_properties(domain):
#    
#    sigma_true=Scalar(sigma_background, Function(domain))
#    Mn_true=Scalar(Mn_background, Function(domain))
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
pde_tol=1e-8
w1=1.e1
usel1=False #or True
epsl1=1e-15
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
core = ['Core']
#   tag(s) for face elements (excluding top surface)
faces = ['faces']
restartfile = 'restart'

true_properties=None
