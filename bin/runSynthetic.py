#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
import argparse
sys.path.append(os.getcwd())
import numpy as np
from fingal import readElectrodeLocations, readSurveyData, makeTagField, setupERTPDE
from esys.finley import ReadMesh

from esys.weipa import saveVTK, saveSilo
from esys.escript.pdetools import Locator, MaskFromTag

parser = argparse.ArgumentParser(description='creates a synthetic survey', epilog="l.gross@uq.edu.au, version Jan 2021")
parser.add_argument('--noise', '-n',  dest='noise', default=0., metavar='NOISE', type=float, help="%% of noise to be added. (default is 0) ")
parser.add_argument('--fullwaver', '-f', dest='fullwaver',  action='store_true', default=False, help='creates a fullwaver-style survey.')
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--silo', '-s',  dest='silofile', metavar='SILO', help="silo file for saving mesh file for visualization (no extension) (output if set).")
parser.add_argument('--plotA', '-A',  dest='plotA', metavar='PLOTA', type=int, default=None, help="electrode A for output (SILO must be set)")
parser.add_argument('--plotB', '-B',  dest='plotB', metavar='PLOTB', type=int, default=None, help="electrode B for output (SILO must be set)")
args = parser.parse_args()


if getMPIRankWorld() == 0: print("** This creates a synthetic survey data set**")

    

config = importlib.import_module(args.config)

    
if 'R' in config.datacolumns:
    usePotentials=True
else:
    usePotentials=False
    
if 'E' in config.datacolumns:
    useFields=True
else:
    useFields=False
    
if any( [ s.find("ERR") >=0 for s in config.datacolumns ]):
    assert args.noise, 'Noise level must be positive.'
    addError=True
else:
    addError= args.noise>0
    
if getMPIRankWorld() == 0: 
    print("configuration "+args.config+" imported.")
    print(f"Data columns to be generated: {config.datacolumns}")
    if usePotentials: print("Potential based data are generated.")
    if useFields: print("Electric field based data are generated.")
    if addError: print(f"Error of {args.noise}% is added to data.")    

domain=ReadMesh(config.meshfile)
if getMPIRankWorld() == 0: print("mesh read from "+config.meshfile)

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
if getMPIRankWorld() == 0: print("%s electrode locations read from %s."%(len(elocations), config.stationfile))
      
survey=readSurveyData(config.schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=[], 
                      dipoleInjections=config.dipoleInjections, 
                      dipoleMeasurements=config.dipoleMeasurements, 
                      delimiter=config.datadelimiter, 
                      commend='#', printInfo=True)

# set the true sigma and gamma:
assert config.true_properties, f"True properties must be defined. See true_properties in {args.config}.py"
sigma_true, gamma_true = config.true_properties(domain) 

txt1=str(sigma_true).replace("\n",';')
txt2=str(gamma_true).replace("\n",';')
if getMPIRankWorld() == 0: 
    print(f"True conductivity sigma_true = {txt1}.")
    print(f"True modifies chargeability gamma_true = {txt2}.")


# set locators to extract predictions:
station_locations=[]
for s in survey.getStationNumeration():
   station_locations.append(survey.getStationLocation(s))

nodelocators=Locator(Solution(domain), station_locations)
elementlocators=Locator(ReducedFunction(domain), station_locations)

if getMPIRankWorld() == 0: print( str(len(station_locations))+ " station locators calculated.")


# PDE:
pde=setupERTPDE(domain)
x=pde.getDomain().getX()[0]
y=pde.getDomain().getX()[1]
z=pde.getDomain().getX()[2]
q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z))
pde.setValue(q=q)


primary_field_field={}
primary_field_solution={}
primary_potential={}
primary_field={}

SIGMA0=config.sigma0
ETA0=config.eta0

if getMPIRankWorld() == 0: print("primary field:")
pde.setValue(A=SIGMA0*kronecker(3))   
for ip in survey.getListOfInjectionStations():
    s=Scalar(0.,DiracDeltaFunctions(pde.getDomain()))
    if config.stationsFMT is None:
        s.setTaggedValue(ip,1.)
    else:    
        s.setTaggedValue(config.stationsFMT%ip,1.)
    pde.setValue(y_dirac=s)
    primary_field_solution[ip]=pde.getSolution()
    primary_field_field[ip]=-grad(primary_field_solution[ip], ReducedFunction(domain))
    primary_potential[ip]=nodelocators(primary_field_solution[ip])
    primary_field[ip]=elementlocators(primary_field_field[ip]) 
    txt1=str(primary_field_solution[ip])
    if getMPIRankWorld() == 0:  print("\t%s : %s "%(ip,txt1))
if getMPIRankWorld() == 0: print(str(len(primary_field_field))+" primary fields calculated.")

def getSecondaryPotentials(sigma_s):
    secondary_field_solution={}
    secondary_field_field={}
    secondary_potential={}
    secondary_field={}
    pde.setValue(A=sigma_s*kronecker(3), y_dirac=Data())
    
    for ip in survey.getListOfInjectionStations():
    
        pde.setValue(X=(SIGMA0-sigma_s)*grad(primary_field_solution[ip]))
        secondary_field_solution[ip]=pde.getSolution()
        secondary_field_field[ip]=-grad(secondary_field_solution[ip], ReducedFunction(domain))
        secondary_potential[ip]=nodelocators(secondary_field_solution[ip])
        secondary_field[ip]=elementlocators(secondary_field_field[ip])         
        txt1=str(secondary_field_solution[ip])
        if getMPIRankWorld() == 0:  print("\t%s : %s "%(ip,txt1))
    return secondary_potential, secondary_field, secondary_field_solution, secondary_field_field

if getMPIRankWorld() == 0: print("secondary  potentials:")
secondary_potential, secondary_field, secondary_field_solution, secondary_field_field = getSecondaryPotentials(sigma_true)
if getMPIRankWorld() == 0: print(str(len(secondary_field_solution))+" secondary potentials calculated.")

if getMPIRankWorld() == 0: print("secondary  potentials chargeability:")
secondary_potential_hat, secondary_field_hat, secondary_field_solution_hat, secondary_field_field_hat = getSecondaryPotentials(sigma_true*1/(1+gamma_true)/(1-ETA0))
if getMPIRankWorld() == 0: print(str(len(secondary_field_field_hat))+" secondary potentials chargeability calculated.")

dVp=survey.makeResistencePrediction(values=primary_potential)
dVs1=survey.makeResistencePrediction(values=secondary_potential)
dVs2=survey.makeResistencePrediction(values=secondary_potential_hat)
dEp=survey.makeResistencePrediction(values=primary_field)
dEs1=survey.makeResistencePrediction(values=secondary_field)
dEs2=survey.makeResistencePrediction(values=secondary_field_hat)

#  now the data file can be created:
if getMPIRankWorld() == 0:

    if args.noise>0:
        std=args.noise/100.
        print("%g %% noise added."%args.noise)
    else:
        std=None
        print("no noise added.")
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
            dVpm=dVp[s+(M,)]-dVp[s+(N,)]
            dVs1m=dVs1[s+(M,)]-dVs1[s+(N,)]
            dVs2m=dVs2[s+(M,)]-dVs2[s+(N,)]
            dEp=dEp[s+(M,)]-dEp[s+(N,)]
            dEp1=dEp1[s+(M,)]-dEp1[s+(N,)]
            dEp2=dEp2[s+(M,)]-dEp2[s+(N,)]
        else:    
            M=m[0]
            if config.usesStationCoordinates:
                out+=FMTX%survey.getStationLocation(M)
            else:
                out+=FMTI%M
            dVpm=dVp[s+(M,)]
            dVs1m=dVs1[s+(M,)]
            dVs2m=dVs2[s+(M,)]
            dEpm=dEp[s+(M,)]
            dEs1m=dEs1[s+(M,)]
            dEs2m=dEs2[s+(M,)]
            
        R=(dVs1m+dVpm)
        Rhat=(dVs2m+dVpm)/(1-ETA0)
        E=(dEs1m+dEpm)
        Ehat=(dEs2m+dEpm)/(1-ETA0)
        EI=sqrt(E[0]**2+E[1]**2+E[2]**2)
        if useFields:
            GAMMA=inner(Ehat-E, E)/EI**2           
        else:
            GAMMA=(R-Rhat)/R
        ETA=GAMMA/(GAMMA+1)
        for o in config.datacolumns:
            if std:
                pert=np.random.normal(0.0, scale=std)
            else:
                pert=0
            if o is 'R': # resistance [V/A]
                out+=FMTG%(R*(1+pert))
            elif o is 'E': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(EI*(1+pert))
            elif o is 'E0': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E[0]*(1+pert))
            elif o is 'E1': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E[1]*(1+pert))
            elif o is 'E2': #  electric field intensity per charge current [V/(Am)]
                out+=FMTG%(E[2]*(1+pert))
            elif o is 'ETA': # chargeability  [1]
                out+=FMTG%(ETA*(1+pert))
            elif o is 'GAMMA': # modified chargeability (=ETA/(1-ETA)) [1]
                out+=FMTG%(GAMMA*(1+pert))
            elif o is 'ERR_R': # resistance [V/A]
                ERR_R=R*std
                out+=FMTG%ERR_R
            elif o is 'ERR_E': # electric field intensity per charge current [V/(Am)] 
                ERR_E=E*std
                out+=FMTG%ERR_E
            elif o is 'ERR_ETA': # chargeability  [1]
                ERR_ETA=ETA*std
                out+=FMTG%ERR_ETA
            elif o is 'ERR_GAMMA': #  modified chargeability (=ETA/(1-ETA)) [1]
                ERR_GAMMA=GAMMA*std
                out+=FMTG%ERR_GAMMA
            elif o is 'RELERR_R': # resistance [1]
                RELERR_R=std
                out+=FMTG%RELERR_R                
            elif o is 'RELERR_E': # electric field intensity per charge current [1]
                ERR_GAMMA=std
                out+=FMTG%ERR_GAMMA
            elif o is 'RELERR_ETA': # chargeability  [1]
                RELERR_ETA=std
                out+=FMTG%RELERR_ETA
            elif o is 'RELERR_GAMMA': #  modified chargeability (=ETA/(1-ETA)) [1]
                RELERR_GAMMA=std
                out+=FMTG%RELERR_GAMMA        
        f.write(out[:-3]+"\n")
        n+=1
    f.close()
    print(n," measurement written to file "+config.datafile)   

if args.silofile is not None:
    sigma_true.expand()
    gamma_true.expand()
    if args.plotA is not None:
        A=int(args.plotA)
    else:
        A=t[0]
    if args.plotB is not None:
        B=int(args.plotB)
    else:
        B=t[1]
    Es=secondary_field_field[A]-secondary_field_field[B]
    Es_hat=(secondary_field_field_hat[A]-secondary_ffield_field_hat[B])/(1-ETA0)
    E=primary_field_field[A]-primary_field_field[B]+Es
    gamma=inner(Es_hat, E)/inner(E, E)

    saveSilo(args.silofile,tag=makeTagField(ReducedFunction(domain)), sigma_true=sigma_true, gamma_true=gamma_true, gamma=gamma, Es=Es, Es_hat=EsHAT,  E=E)
    print(args.silofile+".silo with tags has been generated for injection [%d, %d]"%(A,B))
