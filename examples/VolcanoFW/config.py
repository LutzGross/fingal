meshfile='mesh_volcano.fly'
region_fixed = None
L_stations=0.3
adjustStationLocationsToElementCenter=True
usesStationCoordinates=False
stationdelimiter=' '
stationfile='stations.loc'
stationsFMT=None
datacolumns=['E', 'GAMMA']
datadelimiter=','
datafile='Data.csv'
schedulefile='FWSchedule.csv' # for synthetic cases
dipoleMeasurements=False
dipoleInjections=True
weightLogDefect=0.5
clip_property_function=10

restartfile='sigmarestart'
tolerance=1e-4
imax=500
truncation=20
restart=60
interpolation_order=3
core=[ "R1", "R2"]

region_fixed = None

w1=0
w0=1e-15
alpha0=1.
alpha1=(570)**2/2.


w1gamma=1e-8
w0gamma=0
alpha0gamma=alpha0
alpha1gamma=alpha1

outfile='sigma'

sigma0=0.001
eta0=0.005
sigma_background=sigma0
eta_background=eta0



def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma_true=Scalar(sigma_background, Function(domain))
    sigma_true.setTaggedValue("R1", 0.01)
    sigma_true.setTaggedValue("R2", 0.1)
    
    gamma_background=eta_background/(1-eta_background) 
    gamma_true=Scalar(gamma_background, Function(domain))
    gamma_true.setTaggedValue("R1", 0.02)
    gamma_true.setTaggedValue("R2", 0.2)
   
    return sigma_true, gamma_true
