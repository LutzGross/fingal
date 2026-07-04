#
# this is a fingal configuration file
#
created = 'today'
project='Synthetic Temperature'
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
datacolumns = ['R', 'ERR_R', 'R2', 'ERR_R2', 'ETA']
datacolumns = ['R', 'R2', 'ETA']
dipoleInjections = True
dipoleMeasurements = True
datadelimiter = ','
usesStationCoordinates = False
schedulefile = 'schedule.csv'

surfacetemperature_file = "surface_temperature.csv"
surfacetemperature_undef = -9999
surfacetemperature_skiprows = 1

surfaceflux_file = "surface_flux.csv"
surfaceflux_undef = 999
surfaceflux_file = None

temperature_longrange = 15.5
#
#  This section of the file defines the inversion
#
#def true_properties(domain):
#    from esys.escript import Scalar, Function
#    sigma_true=Scalar(sigma_ref, Function(domain))
#    Mn_true=Scalar(Mn_ref, Function(domain))
#    return sigma_true, Mn_true
#
pde_tol = 1e-10
use_log_misfit_DC=False
data_rtol = 1e-4
use_log_misfit_IP = False
regularization_weighting_DC_misfit = 1
clip_property_function = 10
#
#  Inversion:
#
m_tolerance=1.e-3
g_tolerance=1.e-4
interpolation_order=3
imax=400
truncation=20
restart=60
#
# Output handeling:
#
outfile='sigma'
restartfile = 'restart'

thermal_conductivity = 2
background_heat_production = 0
background_heat_production_core = 0
background_heat_surface_flux_bottom = 0
background_heat_surface_flux_top = 0

# regularization:

model = "TEMP"

useTemperatureExtrapolation = True

regularization_w1DC = 0.005
regularization_w1DC = 0.01
#regularization_w1DC = 0.1
regularization_w1DC = 0.0001
surface_values_variation_factor = 1000