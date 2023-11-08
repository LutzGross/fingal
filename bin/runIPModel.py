#!/usr/bin/python3
#from __future__ import print_function
from esys.escript import *
import importlib, sys, os
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "lib"))
sys.path.append(os.getcwd())
from fingal import IPModel,  readElectrodeLocations, readSurveyData
from esys.finley import ReadMesh, ReadGmsh
from esys.escript import unitsSI as U
import numpy as np
from esys.weipa import saveVTK, saveSilo
import argparse
from specsim3d import spectral_random_field
from esys.escript.pdetools import Locator


parser = argparse.ArgumentParser(description='creates a synthetic ERT/IP survey data set.', epilog="l.gross@uq.edu.au, version 9/8/2023")
parser.add_argument(dest='config', metavar='configfile', type=str, help='python setting configuration')
parser.add_argument('--noise', '-n',  dest='noise', default=0., metavar='NOISE', type=float, help="%% of noise to be added. (default is 0) ")
parser.add_argument('--silo', '-s',  dest='silofile', metavar='SILO', help="silo file for saving mesh file for visualization (no extension) (output if set).")
parser.add_argument('--plot', '-p',  dest='plot', action='store_true', default=False, help="create plot of created properties.")

args = parser.parse_args()


print("** This is a ERT/IP synthetic dataset for true property distribution from config. **")


config = importlib.import_module(args.config)
print("configuration "+args.config+" imported.")


domain=ReadMesh(config.meshfile)
print("mesh read from "+config.meshfile)

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))
      
survey=readSurveyData(config.schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=[], dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements, delimiter=config.datadelimiter)

     
print("%s observations read from %s."%(survey.getNumObservations(), config.schedulefile))
if config.dipoleInjections:
    if config.dipoleMeasurements:
        print("Dipole - Dipole survey")
    else:
        print("Dipole - Pole survey")
else:
    if config.dipoleMeasurements:
        print("Pole - Dipole survey")
    else:
        print("Pole - Pole survey")   

sigma0=config.sigma0
sigma_background=config.sigma_background
    
gamma_background=config.eta_background/(1-config.eta_background)
gamma0=(config.eta0/(1-config.eta0))
    
srclocators=Locator(ContinuousFunction(domain),  [  survey.getStationLocation(ip) for i, ip in enumerate(survey.getListOfInjectionStations()) ])    


model = IPModel(domain, survey=survey, locations=elocations, 
                field_resolution=config.dx, 
                field_origin=config.coreOrigin,
                sigma_background=config.sigma_background, 
                gamma_background=config.eta_background/(1-config.eta_background), 
                padding_tags=config.tagsPadding,
                stationsFMT=config.stationsFMT)



origin=config.coreOrigin
if config.SurveyDim ==2:
    dims=(config.Nx, config.Nz)
else:
    dims=(config.Nx, config.Ny, config.Nz)
spec = spectral_random_field(domainsize = dims, covmod = config.covarianceModel)
    
field_sigma = spec.simnew()
field_gamma = spec.simnew()

np.save(config.field_sigma_file, field_sigma)
np.save(config.field_gamma_file, field_gamma)
print("true fields saved to %s and %s"%(config.field_sigma_file, config.field_gamma_file))




sigma_true=sigma0*np.exp(field_sigma)
gamma_true=gamma0*np.exp(field_gamma)
print("sigma0 = %s"%(str(sigma0)))
print("true sigma min, max, mean= %s, %s, %s"%(np.amin(sigma_true), np.amax(sigma_true), np.mean(sigma_true)))
print("gamma0 = %s"%(str(gamma0)))
print("true gamma min, max, mean= %s, %s, %s"%(np.amin(gamma_true), np.amax(gamma_true), np.mean(gamma_true)))
print("true eta min, max, mean= %s, %s, %s"%(np.amin(gamma_true/(1+gamma_true)), np.amax(gamma_true/(1+gamma_true)), np.mean(gamma_true/(1+gamma_true))))

data=model.runSurvey(model.getAllInjections(), sigma_true, gamma_true)

if args.plot:
        import matplotlib.pylab as plt
        plt.figure()
        plt.imshow(sigma_true)
        plt.colorbar()
        plt.savefig('sigma_true.png')
        plt.clf()
        plt.imshow(gamma_true)
        plt.colorbar()
        plt.savefig('gamma_true.png')
        plt.close()
        print("plot files created.")
        
if args.noise>0:
    std=args.noise/100.
    pert=np.random.normal(0.0, scale=std, size=model.use.shape)
    data[model.use]*=(1+pert[model.use])    
    print("%g %% noise added."%args.noise)
else:
    print("no noise added.")




if survey.hasDipoleInjections():
    if config.usesStationCoordinates:
        FMT=("%10.15e"+config.datadelimiter)*12+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"
    else:
        FMT=("%d"+config.datadelimiter)*4+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"            
else:
    if config.usesStationCoordinates:
        FMT=("%10.15e"+config.datadelimiter)*9+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"
    else:
        FMT=("%d"+config.datadelimiter)*3+"%10.15e"+config.datadelimiter+"%10.15e"+"\n"            
f=open(config.datafile,'w')
n=0
for k, j in enumerate(model.getAllInjections()):
    for s in survey.getObservations(model.getInjection(j)):
        #print(j,s,k,model.dataDCMaps[j][s], model.use[model.dataDCMaps[j][s],k])
        if model.use[model.dataDCMaps[j][s],k]:
            V=data[model.dataDCMaps[j][s],k]
            dV=data[model.dataIPMaps[j][s],k]
            eta=dV/(V+dV)
            ips=model.getInjection(j)
            if survey.hasDipoleInjections():
                if config.usesStationCoordinates:
                    xA=survey.getStationLocation(ips[0])
                    xB=survey.getStationLocation(ips[1])
                    xM=survey.getStationLocation(s[0])
                    xN=survey.getStationLocation(s[1])
                    f.write(FMT%(xA[0], xA[1], xA[2], xB[0], xB[1], xB[2], xM[0], xM[1], xM[2], xN[0], xN[1], xN[2], V, eta ))
                else:
                    f.write(FMT%(ips[0], ips[1], s[0], s[1], V, eta ))

            else:
                if config.usesStationCoordinates:
                    xA=survey.getStationLocation(ips)
                    xM=survey.getStationLocation(s[0])
                    xN=survey.getStationLocation(s[1])
                    f.write(FMT%(xA[0], xA[1], xA[2], xM[0], xM[1], xM[2], xN[0], xN[1], xN[2], V, eta ))
                else:
                    f.write(FMT%(ips, s[0], s[1], V, eta ))
            n+=1
f.close()

print(str(n)+" measurement written to file "+config.datafile)

if args.silofile is not None:
    model.sigma.expand() 
    model.gamma.expand()
    saveSilo(args.silofile, sigma=model.sigma, eta=model.gamma/(1+model.gamma), tag=makeTagMap(ReducedFunction(domain)))
    print("Material properties have been written to "+args.silofile+".silo")
