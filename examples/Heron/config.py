#
# this is a fingal configuration file
#
created = 'today'
project='Heron Line'
# 
#   name of the mesh file
# 
meshfile = 'mesh.fly'
faces_tags = ['faces']
surface_tags = ['surface']
core_tags = ['core']
padding_tags = ['padding']
#
#  file of the location of stations/electrodes. 
#
stationfile = 'stations.csv'
stationdelimiter = ','
stationsFMT = 's%s'
#
# this defines the data file: 
# 
datafile = 'data.csv'
datacolumns = ['R']
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = None
data_rtol = 0.001
#
#  This section of the file defines the inversion
#
sigma0_ref=0.002
Mn_ref=0.01*sigma0_ref
true_properties=None
#
#
#  Inversion:
#
fixed_region_tags=[]
fix_top = False
#weighting_misfit_ERT=0.5
use_robin_condition_in_model = False
clip_property_function=10
m_tolerance=1.e-4
g_tolerance=1.e-6
interpolation_order=3
imax=400
truncation=20
restart=600
pde_tol=1e-10
regularization_w1DC=1e-5
#regularization_w1DC=1e-4
use_L1Norm=False
epsilon_L1Norm=0.01
use_log_misfit_DC = False
regularization_order = 'H2_0' # in ['H1', 'H2', 'H1_0', 'H2_0', 'Gauss', DGauss']

regularization_length_scale = 3
regularization_penalty_factor = 10
# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


