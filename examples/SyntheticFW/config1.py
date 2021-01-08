#     gmsh -3 -o test1.msh test1.geo
# check electrode config!
#    python mkSurvey.py 
#    convertERTGmshToFly.py test1.msh config1
#    runFWsynthetic.py -n 5 config1
# 
# copy test1.fly, test1_data.csv, test1.loc  



meshfile="test1.fly"
restartfile='restart1'
#========================
L_stations=0.01
adjustStationLocationsToElementCenter=True # if set the location of observation stations is moved to the nearest element center. L_stations is then almost ignored
                                           # otherwise values are averaged over all quadrature points in a circle of diameter L_stations around a station are used
# file of the list of measuring electrode locations in the form
# A, x, y, z
# where A (int) is the electrode id and x,y,z (floats) is its location.
# 
usesStationCoordinates=False
stationdelimiter=','
stationfile='test1.loc'
stationsFMT=None          # 

# file of the data in row format
# A, B, S

datadelimiter=','
datafile='test1_data.csv'
schedulefile='test1_schedule.csv' # for synthetic cases 
dipoleInjections=True
dipoleMeasurements=False


sigma_background=0.001
sigma_true={"OuterBox": sigma_background, "InnerBox" : sigma_background,    "Anomaly1" : sigma_background*10,   "Anomaly2" : sigma_background*10, "Anomaly3": sigma_background} # for synthetic case
gamma_background=0.01
gamma_true  ={"OuterBox": gamma_background, "InnerBox" : gamma_background,     "Anomaly1" : gamma_background, "Anomaly2" : gamma_background*20 , "Anomaly3": gamma_background*20} # for synthetic case
sigma0=sigma_background
gamma0=gamma_background
region_fixed = None
test_region = [ "InnerBox", "Anomaly1", "Anomaly2", "Anomaly3" ]

useLogDefect=True
outfile='sigma'
sensitivityfile='sens'
tolerance=1e-3
interpolation_order=1

alpha0=1
alpha1=(300.)**2

w0=5.e-14
w0=1.e-18
w1=0.


alpha0gamma=1
alpha1gamma=(300.)**2
w0gamma=1.e-1
w1gamma=0.

 




