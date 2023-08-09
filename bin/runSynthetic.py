#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
import argparse
sys.path.append(os.getcwd())
import numpy as np
from fingal import readElectrodeLocations, readSurveyData, setupERTPDE, getInjectionPotentials
from esys.finley import ReadMesh

from esys.weipa import saveVTK, saveSilo
from esys.escript.pdetools import Locator, MaskFromTag

parser = argparse.ArgumentParser(description='creates a synthetic ERT/IP survey data set', epilog="l.gross@uq.edu.au, version 9/8/2023")
parser.add_argument('--noise', '-n',  dest='noise', default=0., metavar='NOISE', type=int, help="%% of noise to be added. (default is 0) ")
parser.add_argument('--fullwaver', '-f', dest='fullwaver',  action='store_true', default=False, help='creates a fullwaver-style survey.')
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--silo', '-s',  dest='silofile', metavar='SILO', help="silo file for saving mesh file and property distributions for visualization (no extension) (output if set).")
args = parser.parse_args()

print("** This creates a synthetic survey data set from properties set by config.true_properties **")
config = importlib.import_module(args.config)

    
if 'R' in config.datacolumns:
    usePotentials=True
else:
    usePotentials=False
    
if 'E' in config.datacolumns:
    useFields=True
else:
    useFields=False

if args.noise > 0:
    addNoise=True
    std=args.noise/100.
else:
    std=5/100.
    addNoise = False

print("configuration "+args.config+" imported.")
print(f"Data columns to be generated: {config.datacolumns}")
if usePotentials: print("Potential based data are generated.")

if useFields: print("Electric field based data are generated.")
if addNoise: print(f"Noise of {args.noise}% is added to data.")
print(f"Error of {std*100}% is assumed.")


domain=ReadMesh(config.meshfile)
print("mesh read from "+config.meshfile)

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
if getMPIRankWorld() == 0: print("%s electrode locations read from %s."%(len(elocations), config.stationfile))
      
survey=readSurveyData(config.schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=[], 
                      dipoleInjections=config.dipoleInjections, 
                      dipoleMeasurements=config.dipoleMeasurements, 
                      delimiter=config.datadelimiter, 
                      commend='#', printInfo=True)

# set the true sigma and gamma:
assert config.true_properties , f"True properties must be defined. See true_properties in {args.config}.py"
sigma_0_true, Mn_true = config.true_properties(domain) 

#sigma_0_true.expand()
#Mn_true.expand()

#saveSilo("TT", sigma=sigma_0_true, Mn=Mn_true, tags=makeTagMap(ReducedFunction(domain)))

txt1=str(sigma_0_true).replace("\n",';')
txt2=str(Mn_true).replace("\n",';')

print(f"True conductivity sigma_0_true = {txt1}.")
print(f"True normalised chargeability Mn_true = {txt2}.")

# set locators to extract predictions:
station_locations=[]
for s in survey.getStationNumeration():
   station_locations.append(survey.getStationLocation(s))

nodelocators=Locator(Solution(domain), station_locations)
elementlocators=Locator(ReducedFunction(domain), station_locations)

print( str(len(station_locations))+ " station locators calculated.")

SIGMA_S=config.sigma_background


print("... injection field:")
print("sigma_S= ",SIGMA_S)
injection_potential=getInjectionPotentials(domain, SIGMA_S, survey, stationsFMT=config.stationsFMT)

injection_field={}
injection_potential_at_station={}
injection_field_at_station={}


for ip in survey.getListOfInjectionStations():
    injection_field[ip]=-grad(injection_potential[ip], ReducedFunction(domain))
    injection_potential_at_station[ip]=nodelocators(injection_potential[ip])
    injection_field_at_station[ip]=elementlocators(injection_field[ip]) 
    txt1=str(injection_potential[ip])
    print("\t%s : %s "%(ip,txt1))
print(str(len(injection_field))+" injection fields calculated.")


# PDE:
pde=setupERTPDE(domain)
x=pde.getDomain().getX()[0]
y=pde.getDomain().getX()[1]
z=pde.getDomain().getX()[2]
q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z))
pde.setValue(q=q)

# -div(sigma_S grad(V_i))=S
# V_S = injection potential
# secondary potential -div(sigma_0 grad(V_2) = -div(sigma_S-sigma_0 grad(V_i))
#  DC potenential V_0 = V_S + V_2
print(".... DC potential : (increment from injection)")
print("sigma_0= ",sigma_0_true)
potential_0 = {}
field_0 = {}
potential_0_at_stations = {}
field_0_at_stations = {}
pde.setValue(A=sigma_0_true * kronecker(3), y_dirac=Data())

for ip in survey.getListOfInjectionStations():
    pde.setValue(X=(SIGMA_S-sigma_0_true) * grad(injection_potential[ip]))
    potential_0[ip] = pde.getSolution()
    field_0[ip] = -grad(potential_0[ip], ReducedFunction(domain))
    potential_0_at_stations[ip] = nodelocators(potential_0[ip])
    field_0_at_stations[ip] = elementlocators(field_0[ip])
    txt1 = str(potential_0[ip])
    print("\t%s : %s " % (ip, txt1))
print(str(len(potential_0))+" secondary potentials calculated.")


# secondary potential -div(sigma_oo grad(V_over) = -div(Mn grad(V_DC))
#  ETA= V_over/V_DC
print(".... secondary potentials:")
sigma_oo_true=Mn_true+sigma_0_true
print("sigma_oo= ",sigma_oo_true)
secondary_potential = {}
secondary_field = {}
secondary_potential_at_stations = {}
secondary_field_at_stations = {}
pde.setValue(A=sigma_oo_true * kronecker(3), y_dirac=Data())

for ip in survey.getListOfInjectionStations():
    pde.setValue(X=Mn_true * grad(potential_0[ip]+injection_potential[ip]))
    secondary_potential[ip] = pde.getSolution()
    secondary_field[ip] = -grad(secondary_potential[ip], ReducedFunction(domain))
    secondary_potential_at_stations[ip] = nodelocators(secondary_potential[ip])
    secondary_field_at_stations[ip] = elementlocators(secondary_field[ip])
    txt1 = str(secondary_potential[ip])
    print("\t%s : %s " % (ip, txt1))
print(str(len(secondary_potential))+" M potentials calculated.")

dV_i=survey.makeResistencePrediction(values=injection_potential_at_station)
dV_0=survey.makeResistencePrediction(values=potential_0_at_stations)
dV_2=survey.makeResistencePrediction(values=secondary_potential_at_stations)
dE_i=survey.makeResistencePrediction(values=injection_field_at_station)
dE_0=survey.makeResistencePrediction(values=field_0_at_stations)
dE_2=survey.makeResistencePrediction(values=secondary_field_at_stations)

#  now the data file can be created:
if getMPIRankWorld() == 0:
    n=0
    f=open(config.datafile,'w')
    
    FMTX="%g"+config.datadelimiter+" %g"+config.datadelimiter+" %g"+config.datadelimiter
    FMTG="%g"+config.datadelimiter+" "    
    FMTI="%d"+config.datadelimiter+" "
    for t in survey.tokenIterator():

        out=""
        if survey.hasDipoleInjections():
            A=t[0]
            B=t[1]
            if config.usesStationCoordinates:
                out=FMTX%survey.getStationLocation(A)+FMTX%survey.getStationLocation(B)
            else:
                out=FMTI%A + FMTI%B
            m=t[2:]
            s=t[:2]
        else:    
            A=t[0]
            s=t[:1]
            m=t[1:]
            if config.usesStationCoordinates:
                out=FMTX%survey.getStationLocation(A)
            else:
                out=FMTI%A

        if survey.hasDipoleMeasurements():
            M=m[0]
            N=m[1]
            if config.usesStationCoordinates:
                out+=FMTX%survey.getStationLocation(M)+FMTX%survey.getStationLocation(N)
            else:
                out+=FMTI%M + FMTI%N
            #dV_i_M=dV_i[s+(M,)]-dV_i[s+(N,)]
            #dV_0_M=dV_0[s+(M,)]-dV_0[s+(N,)]
            #dV_2_M=dV_2[s+(M,)]-dV_2[s+(N,)]
            #dE_i_M=dE_i[s+(M,)]-dE_i[s+(N,)]
            #dE_0_M=dE_0[s+(M,)]-dE_0[s+(N,)]
            #dE_2_M=dE_2[s+(M,)]-dE_2[s+(N,)]
        else:    
            M=m[0]
            if config.usesStationCoordinates:
                out+=FMTX%survey.getStationLocation(M)
            else:
                out+=FMTI%M
            #dV_i_M=dV_i[s+(M,)]
            #dV_0_M=dV_0[s+(M,)]
            #dV_2_M=dV_2[s+(M,)]
            #dE_i_M=dE_i[s+(M,)]
            #dE_0_M=dE_0[s+(M,)]
            #dE_2_M=dE_2[s+(M,)]
        V_0 =dV_0[t]+dV_i[t]
        V_oo=dV_0[t]+dV_i[t]-dV_2[t]
        E_0=dE_0[t]+dE_i[t]
        E_oo=dE_0[t]+dE_i[t]-dE_2[t]
        EI=sqrt(E_0[0]**2+E_0[1]**2+E_0[2]**2)
        ETA = dV_2[t] / V_0
        if useFields:
            raise ValueError("CHECK GAMMA")
            GAMMA=inner(E_oo-E_0, E_0)/EI**2           
        else:
            GAMMA=ETA/(1-ETA)
        
        for o in config.datacolumns:
            if addNoise:
                pert=np.random.normal(0.0, scale=std)
            else:
                pert=0
            if o == 'R': # resistance [V/A]
                out+=FMTG%(V_0*(1+pert))
            elif o == 'E': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(EI*(1+pert))
            elif o == 'E0': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E_0[0]*(1+pert))
            elif o == 'E1': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E_0[1]*(1+pert))
            elif o == 'E2': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E_0[2]*(1+pert))
            elif o == 'ETA': # chargeability  [1]
                out+=FMTG%(ETA*(1+pert))
            elif o == 'GAMMA': # modified chargeability (=ETA/(1-ETA)) [1]
                raise ValueError("check gamma")
                out+=FMTG%(GAMMA*(1+pert))
            elif o == 'ERR_R': # resistance [V/A]
                print(V_0, std, pert)
                ERR_R=abs(V_0)*std
                out+=FMTG%ERR_R
            elif o == 'ERR_E': # electric field intensity per charge current [V/(Am)] 
                ERR_E=abs(EI)*std
                out+=FMTG%ERR_E
            elif o == 'ERR_ETA': # chargeability  [1]
                ERR_ETA=abs(ETA)*std
                out+=FMTG%ERR_ETA
            elif o == 'ERR_GAMMA': #  modified chargeability (=ETA/(1-ETA)) [1]
                raise ValueError("check gamma")
                ERR_GAMMA=abs(GAMMA)*std
                out+=FMTG%ERR_GAMMA
            elif o == 'RELERR_R': # resistance [1]
                RELERR_R=std
                out+=FMTG%RELERR_R                
            elif o == 'RELERR_E': # electric field intensity per charge current [1]
                RELERR_E=std
                out+=FMTG%RELERR_E
            elif o == 'RELERR_ETA': # chargeability  [1]
                RELERR_ETA=std
                out+=FMTG%RELERR_ETA
            elif o == 'RELERR_GAMMA': #  modified chargeability (=ETA/(1-ETA)) [1]
                RELERR_GAMMA=std
                out+=FMTG%RELERR_GAMMA        
        f.write(out[:-3]+"\n")
        n+=1
    f.close()
    print(n," measurement written to file "+config.datafile)   

if args.silofile is not None:
    sigma_0_true.expand()
    Mn_true.expand()
    saveSilo(args.silofile,tag=makeTagMap(ReducedFunction(domain)), sigma_0_true=sigma_0_true, Mn_true=Mn_true)
    print(args.silofile+".silo with tags has been generated for injection [%d, %d]"%(A,B))
