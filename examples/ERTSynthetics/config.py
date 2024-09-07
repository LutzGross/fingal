#
# this is a fingal configuration file
#
created = 'today'
project='ERTSynthetics'
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
#
#  This section of the file defines the inversion
#
#region_fixed=[]
sigma0_ref=1e-5
Mn_ref=0.01*sigma0_ref
def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma0_true=Scalar(sigma0_ref, Function(domain))
    sigma0_true.setTaggedValue('anomaly', sigma0_ref * 100)
    Mn_true=Scalar(Mn_ref, Function(domain))
    Mn_true.setTaggedValue('anomaly', sigma0_ref * 100 * 0.25)
    return sigma0_true, Mn_true
#
#
#  Inversion:
#
#weighting_misfit_ERT=0.5
#clip_property_function=10
#m_tolerance=1.e-4
#g_tolerance=1.e-4
#interpolation_order=3
#imax=400
#truncation=20
#restart=60
#pde_tol=1e-10
#w1=0.01
#w1=1e-8
#useL1Norm=False
#epsilonL1Norm=1e-4
#useLogMisfit = False

# Output handeling:
#
#outfile='sigma'

#restartfile = 'restart'


