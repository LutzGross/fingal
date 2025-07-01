#
# this is a fingal configuration file
#
created = 'today'
project='IPSynthetics'
# 
#   name of the mesh file
# 
meshfile = 'mesh.fly'
meshfile = 'mesh_synthetic.fly'
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
datacolumns = ['R', 'ETA']
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = 'schedule.csv'
data_rtol = 1e-4
#
#  This section of the file defines the inversion
#
sigma0_ref=0.002
Mn_ref=0.01*sigma0_ref

sigma0_ref=0.1
Mn_ref=0.001*sigma0_ref

def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma0_true=Scalar(sigma0_ref , Function(domain))
    sigma0_true.setTaggedValue('anomaly', sigma0_ref * 10)
    Mn_true=Scalar(Mn_ref, Function(domain))
    Mn_true.setTaggedValue('anomaly',  sigma0_ref * 10 * 0.5 )
    return sigma0_true, Mn_true

#
#
#  Inversion:
#
use_robin_condition_in_model = False
use_L1Norm = False
epsilon_L1Norm = False
fixed_region_tags=[]
fix_top = False
#weighting_misfit_ERT=0.5
clip_property_function=10
m_tolerance=1.e-4
g_tolerance=1.e-6
interpolation_order=3
imax=400
truncation=20
restart=60
pde_tol=1e-10
regularization_w1DC=1e-1
regularization_w1IP=1e-1
regularization_theta = 0.
#regularization_w1DC=1e-4
use_log_misfit_DC = False
use_log_misfit_IP = False

regularization_weighting_DC_misfit =  1
regularization_order = 'H1' # in ['H1', "H1_0", 'H2',  "H2_0"]
if regularization_order == 'H2_0' :
    regularization_length_scale = None # only used for "H2" and "H2_0" regularization
    regularization_w1DC=100
    regularization_w1IP=100
if regularization_order == 'H1_0' :
    regularization_w1DC=0.001
    regularization_w1IP = 0.001
# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


