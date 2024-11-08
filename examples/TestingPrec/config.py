#
# this is a fingal configuration file
#
created = 'today'
project='testing'
# 
#   name of the mesh file
# 
meshfile = 'brick.fly'
faces_tags = ['faces']
surface_tags = []
core_tags = ['core', 'anomaly']
padding_tags = []
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
schedulefile = 'schedule.csv'
data_rtol = 1e-4
#
#  This section of the file defines the inversion
#
sigma0_ref=0.0004
Mn_ref=0.01*sigma0_ref
def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma0_true=Scalar(sigma0_ref , Function(domain))
    sigma0_true.setTaggedValue('anomaly', sigma0_ref * 100)
    Mn_true=Scalar(Mn_ref, Function(domain))
    Mn_true.setTaggedValue('anomaly', Mn_ref * 100 * 0.25)
    return sigma0_true, Mn_true
#
#
#  Inversion:
#
fixed_region_tags=[]
fix_top = False
#weighting_misfit_ERT=0.5
clip_property_function=10
m_tolerance=1.e-4
g_tolerance=1.e-5
interpolation_order=3
imax=400
truncation=20
restart=6000
pde_tol=1e-10

#regularization_w1=1e-4
use_L1Norm=False
epsilon_L1Norm=0.01
use_log_misfit_DC = False
regularization_order = 'H2' # in ['H1', 'H2', 'Gauss', 'PseudoGauss', D-PseudoGauss']
regularization_length_scale = 5000
regularization_w1=1e-1
# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


