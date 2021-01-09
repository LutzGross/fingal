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





primary_field={}
primary_potential={}
primary_field_solution={}
primary_field_elements={}



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
    primary_potential[ip]=pde.getSolution()
    primary_field[ip]=-grad(primary_potential[ip], ReducedFunction(domain))
    primary_field_solution[ip]=nodelocators(primary_potential[ip])
    primary_field_elements[ip]=elementlocators(primary_field[ip]) 
    txt1=str(primary_potential[ip])
    if getMPIRankWorld() == 0:  print("\t%s : %s "%(ip,txt1))
if getMPIRankWorld() == 0: print(str(len(primary_field))+" primary fields calculated.")

def getSecondaryPotentials(sigma_s):
    secondary_field={}
    secondary_potential={}
    secondary_field_solution={}
    secondary_field_elements={}
    pde.setValue(A=sigma_s*kronecker(3), y_dirac=Data())
    
    for ip in survey.getListOfInjectionStations():
    
        pde.setValue(X=(sigma_s-SIGMA0)*grad(primary_potential[ip]))
        secondary_potential[ip]=pde.getSolution()
        secondary_field[ip]=-grad(secondary_potential[ip], ReducedFunction(domain))
        if usePotentials:
            secondary_field_solution[ip]=nodelocators(secondary_potential[ip])
        if useFields:
            secondary_field_elements[ip]=elementlocators(secondary_field[ip])         
        txt1=str(secondary_potential[ip])
        if getMPIRankWorld() == 0:  print("\t%s : %s "%(ip,txt1))
    return secondary_field, secondary_potential, secondary_field_solution, secondary_field_elements

if getMPIRankWorld() == 0: print("secondary  potentials:")
secondary_field, secondary_potential, secondary_field_solution, secondary_field_elements = getSecondaryPotentials(sigma_true)
if getMPIRankWorld() == 0: print(str(len(secondary_field))+" secondary potentials calculated.")

if getMPIRankWorld() == 0: print("secondary  potentials chargeability:")
secondary_field_hat, secondary_potential_hat, secondary_field_solution_hat, secondary_field_elements_hat = getSecondaryPotentials(sigma_true*gamma_true/(1+gamma_true)/ETA0)
if getMPIRankWorld() == 0: print(str(len(secondary_field_hat))+" secondary potentials chargeability calculated.")

1/0
#======================================
if hasattr(config, 'dipoles'):
    dipoles=config.dipoles
else:
    dipoles=False
if dipoles:
    print("Dipoles at stations are used.")






# set the reference conductivity:
sigma_true=Scalar(config.sigma_background,Function(domain))
gamma_true=Scalar(config.gamma_background,Function(domain))
for k in config.sigma_true:
    sigma_true.setTaggedValue(k, config.sigma_true[k])
    gamma_true.setTaggedValue(k, config.gamma_true[k])

if isinstance(config.sigma0, dict):
    sigma0=Scalar(config.sigma_background,Function(domain))
    for k in config.sigma0:
        sigma0.setTaggedValue(k, config.sigma0[k])
else:
    sigma0=config.sigma0

if isinstance(config.gamma0, dict):
    gamma0=Scalar(config.gamma_background,Function(domain))
    for k in config.sigma0:
        gamma0.setTaggedValue(k, config.gamma0[k])
else:
    gamma0=config.gamma0
   
print("sigma 0 = %s"%(str(sigma0)))
print("true sigma = %s"%(str(sigma_true)))
print("gamma 0 = %s"%(str(gamma0)))
print("true gamma = %s"%(str(gamma_true)))




print("secondary potential:")
secondary_field={}
secondary_fieldHAT={}

secondary_potential={}
secondary_potentialHAT={}
secondary_field_loc={}
secondary_fieldHAT_loc={}

pde.setValue(A=sigma_true*kronecker(3), y_dirac=Data())
pdeHAT.setValue(A=sigma_true*kronecker(3), y_dirac=Data())   

for ip in survey.getListOfInjectionStations():
    
    pde.setValue(X=(sigma0-sigma_true)*grad(primary_potential[ip]))
    secondary_potential[ip]=pde.getSolution()
    secondary_field[ip]=-grad(secondary_potential[ip], ReducedFunction(domain))
    secondary_field_loc[ip]=locators(secondary_field[ip])
    
    pdeHAT.setValue(X=sigma_true*gamma_true*grad(primary_potential[ip]+secondary_potential[ip]) )
    secondary_potentialHAT[ip]=pdeHAT.getSolution()
    secondary_fieldHAT[ip]=-grad(secondary_potentialHAT[ip], ReducedFunction(domain))
    secondary_fieldHAT_loc[ip]=locators(secondary_fieldHAT[ip])


    print( "%s : %s : %s"%(ip, secondary_potential[ip], secondary_potentialHAT[ip]))
print(str(len(secondary_field))+" secondary fields calculated.")

dVp=survey.makeResistencePrediction(values=primary_field_loc)
dVs=survey.makeResistencePrediction(values=secondary_field_loc)
dVsHAT=survey.makeResistencePrediction(values=secondary_fieldHAT_loc)

# create a data vector
dataindex={}
data=[]
for A,B, M in dVp:
    xA=survey.getStationLocation(A)
    xB=survey.getStationLocation(B)
    xM=survey.getStationLocation(M)
    dataindex[(A,B,M)]=len(data)

    if dipoles:
        E=dVp[(A,B,M)]+dVs[(A,B,M)]
        dE=dVsHAT[(A,B,M)]-dVs[(A,B,M)]
        data.append((E[0], E[1], 0., dE[0], dE[1], 0.))
    else:
        E= length(dVp[(A,B,M)]+dVs[(A,B,M)])  
        gamma=inner(dVsHAT[(A,B,M)], dVs[(A,B,M)]+dVp[(A,B,M)])/E**2     
        data.append((E, gamma))
# add noise:

data=np.array(data)
if args.noise>0:
    std=args.noise/100.
    for c in range(data.shape[1]):
        pert=np.random.normal(0.0, scale=std, size=data.shape[0])
        data[:,c]*=(1+pert)    
    print("%g %% noise added."%args.noise)
else:
    print("no noise added.")
    

f=open(config.datafile,'w')
if config.usesStationCoordinates:
    if dipoles:
        FMT=("%10.15e"+config.datadelimiter)*9+("%10.15e"+config.datadelimiter)*5+" %10.15e"+"\n"
        for A,B, M in dataindex:
            xA=survey.getStationLocation(A)
            xB=survey.getStationLocation(B)
            xM=survey.getStationLocation(M)        
            f.write(FMT%(xA[0], xA[1], xA[2], xB[0], xB[1], xB[2], xM[0], xM[1], xM[2], data[dataindex[(A,B,M)],0],data[dataindex[(A,B,M)],1], data[dataindex[(A,B,M)],2], data[dataindex[(A,B,M)],3],data[dataindex[(A,B,M)],3],data[dataindex[(A,B,M)],5]))
    else:
        FMT=("%10.15e"+config.datadelimiter)*9+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"
        for A,B, M in dataindex:
            xA=survey.getStationLocation(A)
            xB=survey.getStationLocation(B)
            xM=survey.getStationLocation(M)        
            f.write(FMT%(xA[0], xA[1], xA[2], xB[0], xB[1], xB[2], xM[0], xM[1], xM[2], data[dataindex[(A,B,M)],0],data[dataindex[(A,B,M)],1]  ))
else:
    if dipoles:
        FMT=("%d"+config.datadelimiter)*3+("%10.15e"+config.datadelimiter)*5+"%10.15e"+"\n"
        for A,B, M in dataindex:
            f.write(FMT%(A,B,M, data[dataindex[(A,B,M)],0], data[dataindex[(A,B,M)],1], data[dataindex[(A,B,M)],2], data[dataindex[(A,B,M)],3], data[dataindex[(A,B,M)],4], data[dataindex[(A,B,M)],5]))
    else:
        FMT=("%d"+config.datadelimiter)*3+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"
        for A,B, M in dataindex:
            f.write(FMT%(A,B,M, data[dataindex[(A,B,M)],0],data[dataindex[(A,B,M)],1] ))

print(str( len(dVp))+" measurement written to file "+config.datafile)

if args.silofile is not None:
    sigma_true.expand()
    A0, B0, M = list(dVp.keys())[0]
    if args.plotA is not None:
        A=int(args.plotA)
    else:
        A=A0
    if args.plotB is not None:
        B=int(args.plotB)
    else:
        B=B0
    Es=secondary_field[A]-secondary_field[B]
    EsHAT=secondary_fieldHAT[A]-secondary_fieldHAT[B]
    E=primary_field[A]-primary_field[B]+Es
    
    gamma=inner(EsHAT, E)/inner(E, E)

    saveSilo(args.silofile,tag=makeTagField(ReducedFunction(domain)), sigma_true=sigma_true, gamma=gamma, E_secondary=Es, Ehat_secondary=EsHAT,  E=E)
    print(args.silofile+".silo with tags has been generated for injection [%d, %d]"%(A,B))
