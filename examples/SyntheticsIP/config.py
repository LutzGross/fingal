#
# this is a fingal configuration file
#
created = 'today'
project='IPSynthetics'
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
def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma0_true=Scalar(sigma0_ref , Function(domain))
    sigma0_true.setTaggedValue('anomaly_left', sigma0_ref * 100)
    sigma0_true.setTaggedValue('anomaly_right', sigma0_ref / 100)
    Mn_true=Scalar(Mn_ref, Function(domain))
    Mn_true.setTaggedValue('anomaly_left', Mn_ref * 100 * 0.25)
    Mn_true.setTaggedValue('anomaly_right', Mn_ref * 0.25)
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
g_tolerance=1.e-6
interpolation_order=3
imax=400
truncation=20
restart=60
pde_tol=1e-10
regularization_w1=1e-2
regularization_theta = 0.
#regularization_w1=1e-4
use_log_misfit_DC = False
use_log_misfit_IP = False
regularization_weighting_DC_misfit =  1
regularization_order = 'H1' # in ['H1', "H1_0", 'H2',  "H2_0"]
regularization_length_scale = 3 # only used for "H2" and "H2_0" regularization
# Output handeling:
#
outfile='sigma'

#restartfile = 'restart'


