#
# this is a fingal configuration file
#
created = '08.01.2021 18:49'
#. 
#. This is the project name. This name is used when ever a file name is needed and
#. no specific name was set.
project='synth'
# 
#   name of the mesh file
# 
#. This is mesh file in the fly format. It is assumed that all station/electrodes are nodes of the mesh
#. and all charging electrodes are added as nodal elements with appropriate tagging.
#.  
meshfile = 'synth.fly'
#
#  file of the location of stations/electrodes. 
#
#.  This is a CSV file. Each line gives a identifier (integer) of the station and the x,y,z coordinates (floats).
#.
stationfile = 'stations.csv'
#.
#. This is the delimiter between columns in the station file: 
#.
stationdelimiter = ','
#.
#. This defines how the station identifier is converted into the station tag used 
#. in the fly file for nodal elements. If None the identifier is used as integer tag when
#. referring to the nodal elements typically when defining charges.
#. 
stationsFMT = 's%s'
#.
#. length scale at stations in meter which defines the radius over which prediction are 
#. averaged to compare with data when predictions are not available at stations. This should
#. be in the order of the element size near stations.
#. 
L_stations=1.
#.
#. If True adjustment of station location is allowed to match mesh positions.
#.
adjustStationLocationsToElementCenter=True
#
# this defines the data file: 
# 
#. The data file is a CSV where each line is of the form 
#.
#.     A, B, M, N, R, GAMMA
#.
#. where A,B is the charging dipole, M,N is measuring dipole, R is the resistance (=potential difference across M,N divided by current),
#. and GAMMA is the modified chargeability which is defined as GAMMA=ETA/(1-ETA) where ETA is conventional chargeability. GAMMA is more suitable
#. for inversion as its values are unbounded.
#.
#. If errors are given then the format is
#.
#.     A, B, M, N, R, RELERR_R, GAMMA, RELERR_GAMMA
#.
#. where RELERR_R is the relative error for R and RELERR_GAMMA is the relative error for GAMMA. 
#. For electric field inversion R is to be read as field intensity.
#. 
datafile = 'synth_data.csv'
#.
#. This defines the data columns of the data file.  
#. It is a list of the strings giving the meaning of the column:
#.
#.   * 'R' - resistance [V/A]
#.   * 'E' - electric field intensity per charge current [V/(Am)] 
#.   * 'ETA' - chargeability  [1]
#.   * 'GAMMA' - modified chargeability (=ETA/(1-ETA)) [1]
#.   * 'ERR_R' - resistance [V/A]
#.   * 'ERR_E' - electric field intensity per charge current [V/(Am)] 
#.   * 'ERR_ETA' - chargeability  [1]
#.   * 'ERR_GAMMA' - modified chargeability (=ETA/(1-ETA)) [1]
#.   * 'RELERR_R' - resistance [1]
#.   * 'RELERR_E' - electric field intensity per charge current [1] 
#.   * 'RELERR_ETA' - chargeability  [1]
#.   * 'RELERR_GAMMA' - modified chargeability (=ETA/(1-ETA)) [1]
#.
datacolumns = ['E', 'ETA']
#. 
#. If False the injections are assumed to pole injections and column B is dropped in the data file. 
#.
dipoleInjections = True
#. 
#. If False the measurements are assumed to be pole measurements and column N is dropped in the data file. 
#. This the expected for field intensity surveys.
#.
dipoleMeasurements = False
#.
#. column delimiter in the data file 
#.
datadelimiter = ','
#.
#. if True station coordinates rather then station identifiers are used in the data files (not recommended)
#.
usesStationCoordinates = False
#.
#. This file is used when generating synthetic data. The file contains lines of the for A, B, M, N (or similar)
#. to define charging dipole (A,B) (or pole A) and measuring dipole (M,N)  (or pole M).
#.
schedulefile = 'synth_schedule.csv'
#
#  This section of the file defines the inversion
#
#.
#. reference conductivity (property function p = log(sigma/sigma0))
#.
sigma0=0.001
#.
#. This defines the region where the conductivity and chargeability is not been updated during inversion, for instance 
#. region_fixed=["padding"] blocks updates in the padding region. If region_fixed is not a list (for instance None) 
#. values are fixed at the bottom of the domain. This is required if w0=0 and w1>0.
#.
region_fixed=[]
#.
#. conductivity outside the region of interest where the inversion happens:
#.
sigma_background=sigma0
#.
#. reference chargeability (property function p = log(eta/eta0))
#.   
eta0=0.01
#.
#. chargeability outside the region of interest where the inversion happens:
#.
eta_background=eta0
#.
#. The function defines the true conductivity and modified chargeability. 
#. Obviously this function is used when synthetic data are constructed or for testing and validation.
#. 
def true_properties(domain):
    from esys.escript import Scalar, Function
    sigma_true=Scalar(sigma_background, Function(domain))
    sigma_true.setTaggedValue("InnerBox", sigma_background)
    sigma_true.setTaggedValue("Anomaly1", sigma_background*10)
    sigma_true.setTaggedValue("Anomaly2", sigma_background*10)
    sigma_true.setTaggedValue("Anomaly3", sigma_background)

    gamma_background=eta_background/(1-eta_background)
    gamma_true=Scalar(gamma_background, Function(domain))
    gamma_true.setTaggedValue("InnerBox", gamma_background)
    gamma_true.setTaggedValue("Anomaly1", gamma_background)
    gamma_true.setTaggedValue("Anomaly2", gamma_background*20)
    gamma_true.setTaggedValue("Anomaly3", gamma_background*20)
    
    return sigma_true, gamma_true
#
#  Inversion:
#
#.
#. If set a data misfit is calculated as combination of log scale ( = ( (log(r)-log(R))/RELERR )**2 and
#. quadratic the misfit calculated as ( = ( (r-R)/(R*RELERR) )**2 where R is data, r is prediction and RELERR is relative errorr.
#. weightLogDefect is the weighting for the log term (not supported by all inversions yet.)                                       
#. 
weightLogDefect=0.5
#.
#. values in the property function (in essence this is log(sigma/sigma0)) which are smaller than -clip_property_function
#. and larger then clip_property_function are clipped in order to avoid overflow and underflow in the exponential function.
#.
clip_property_function=10
#.
#. tolerence of the inversion. 
#. Iteration is terminated when both the relative in the property function and the change 
#. in the cost function is less than tolerance.
#.
tolerance=1.e-4
#.
#. Parameters to control Minimizer:
#.
interpolation_order=3
imax=400
truncation=20
restart=60
#.
#.  
#. weighting for the |grad(m)|^2 term
#.
w1=0
#. 
#. weighting for the |m|^2 term
#.
w0=0.1
#.
#. property function p is the smoothed version of m. It is constructed by 
#. solving 
#.     alpha1*p_{,xx} +alpha0*p=m 
#. which applies smmothing with length scale sqrt(alpha1)/alpha0 and a 
#. rescaling by the factor 1/alpha0. Typically sqrt(alpha1) is set to the distance of stations.
#. 
alpha0=1
alpha1=(300.)**2
#.
#. this is the same for the modified chargeability GAMMA
#.
w1gamma=w1
w0gamma=w0
alpha0gamma=alpha0
alpha1gamma=alpha1
#
# Output handeling:
#
#.
#. This is name of the result file. An appropriate extension for the used file format will be added.
#.
outfile='sigma'
#.
#. This defines the core region/region of interest via tags. For instance this region is used to focus output. 
#. If set to to None no selection is applied and the entire domain is used. 
#.
core=["Anomaly1", "Anomaly2", "Anomaly3", "InnerBox"]
#.
#. For an IP inversion the calculation can be split into the conductivity and chargeability part while 
#. the intermediate result is saved to a restart file to allow for a faster configuration of the chargeability
#. inversion. When inversion is running under MPI one restart file per rank is created. This sets the file name
#. for teh restart file:
#.
restartfile = 'restart'

