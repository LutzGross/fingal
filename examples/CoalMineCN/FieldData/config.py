#
# this is a fingal configuration file
#
created = 'today'
project='Coal Mine'
# 
#   name of the mesh file
# 
meshfile = '../mine.fly'
faces_tags = ['Faces']
surface_tags = []
core_tags = ['Seam', 'Goaf', 'Core', 'Base']
padding_tags = ['Padding']
#
#  file of the location of stations/electrodes. 
#
stationfile = 'stations.csv'
stationdelimiter = ','
stationsFMT = 's%s'
#
# this defines the data file: 
# 
datafile = 'data_field.csv'
datacolumns = ['R']
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = '../schedule.csv'
data_rtol = 1e-4
#
#  This section of the file defines the inversion
#
sigma0_ref= 0.025# 0.0188 # 0.0275 # 0.085 # 0.0275 # 0.03 # 0.02753
Mn_ref=0.01*sigma0_ref
true_properties=None
use_robin_condition_in_model = False
#
#
#  Inversion:
#
fixed_region_tags=[]
fix_top = True
#weighting_misfit_ERT=0.5
clip_property_function=10
m_tolerance=1.e-4
g_tolerance=1.e-5
interpolation_order=3
imax=300
truncation=20
restart=6000
pde_tol=1e-10


#regularization_w1DC=1e-4
use_L1Norm=False
epsilon_L1Norm=0.01
use_log_misfit_DC = False


regularization_DC = 'H2' # in ['H1', 'H2', 'H1_0', 'H2_0', 'Gauss', DGauss']

# H2 & H2_0
regularization_w1DC=8e-4
#

# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


