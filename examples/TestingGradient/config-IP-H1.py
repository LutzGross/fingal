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
datacolumns = ['R', "ERR_R", "ETA", "RELERR_ETA"]
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = None
data_rtol = 0.0001
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
pde_tol = 1e-12
fix_top = False
fixed_region_tags=[]
#weighting_misfit_ERT=0.5
clip_property_function=10
regularization_w1=1e-5
#regularization_w1=1e-4
use_L1Norm=False
epsilon_L1Norm=0.01
use_log_misfit_DC = False
regularization_order = 'H1' # in ['H1', 'H2', 'Gauss', 'PseudoGauss', D-PseudoGauss']
regularization_length_scale = 1.3*2



