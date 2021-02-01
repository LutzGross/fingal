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
useLogDefect=True
sigma_background=0.01

restartfile='sigmarestart'
tolerance=1e-4
region_fixed = None
useLogDefect=True

w1=0
w0eta=0
alpha1=(570)**2/2.


alpha0=1.
w0=1e-13

useLogDefectForEta=False
useDifferenceOfFields=True


alpha1eta=95000.
alpha0eta=1.
w1eta=1.e-4
beta=0.1

outfile='sigma'
#eta0=0.001
#sigma0=0.01
#sigma_true={1: sigma_background, 2 : 0.01, 3 : 0.1}
#eta_true={1: eta_background, 2 : 0.001, 3 : 0.2}

#w0=0
#w0eta=0
#alpha1=95000
#alpha0=1.
#w1=5e-8
#useLogDefectForEta=False
#useDifferenceOfFields=True
#alpha1eta=95000.
#alpha0eta=1.
#w1eta=1.e-4 
#beta=0.1
#outfile='sigma'
interpolation_order=1
