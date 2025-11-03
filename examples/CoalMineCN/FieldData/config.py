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
stationfile = '../stations.csv'
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
sigma0_ref= 0.03
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


regularization_DC = 'H2'

CASE = 1

if CASE == 1 :
    regularization_length_scale = None
    regularization_w1DC = 1e-3
    sigma0_ref = 0.03
elif CASE == 2:
    regularization_length_scale = 100
    regularization_w1DC = 1.e-2
    sigma0_ref = 0.03
elif CASE == 3:
    regularization_length_scale = 50
    regularization_w1DC = 1.e-2
    sigma0_ref = 0.03
elif CASE == 4 :
    regularization_length_scale = 5
    regularization_w1DC = 1e-3
    sigma0_ref = 0.03
elif CASE == 5 :
    regularization_length_scale = 1
    regularization_w1DC = 1e-4
    sigma0_ref = 0.03

# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


