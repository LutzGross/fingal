from esys.escript import *
import numpy as np
from math import floor
from scipy.interpolate import RegularGridInterpolator
from .datamapping import mapToDomain
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator

def setupERTPDE(domain, poisson=True):
    """
    used t setup all ERT PDEs
    """
    pde=LinearSinglePDE(domain, isComplex=False)
    pde.setSymmetryOn()
    optionsG=pde.getSolverOptions()
    #optionsG.setSolverMethod(SolverOptions.DIRECT)

    optionsG.setSolverMethod(SolverOptions.PCG)
    optionsG.setTolerance(1e-8)
    if True and hasFeature('trilinos'):
        #print("trilinos solver used.")
        optionsG.setPackage(SolverOptions.TRILINOS)
        optionsG.setPreconditioner(SolverOptions.AMG)
        if poisson:
            optionsG.setTrilinosParameter("problem:type", "Poisson-3D")
        optionsG.setTrilinosParameter("verbosity", "none")
        optionsG.setTrilinosParameter("number of equations", 1)
        #optionsG.setTrilinosParameter("max levels", 3)  # 10 is default 3 seems to be a good number
        #optionsG.setTrilinosParameter("cycle type", "V")
        optionsG.setTrilinosParameter("problem: symmetric", True)
        #optionsG.setTrilinosParameter("smoother: pre or post", "both")
        #optionsG.setTrilinosParameter("Convergence Tolerance", 1e-12)
    return pde


class IPModel(object):
    """
    """
    def __init__(self, domain, survey, locations=[], field_resolution=1., field_origin=(0.,0.,0), sigma_background=0.1, gamma_background=0.0001, padding_tags=[], stationsFMT=None):
        self.domain=domain
        self.survey=survey
        self.locations=locations
        self.stationsFMT=stationsFMT
        self.pde=setupERTPDE(domain)
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.locations=locations
        self.observation_locator=Locator(Solution(domain), [ self.survey.getStationLocation(s) for s in self.survey.getObservationElectrodes()])
        self.source_locator=Locator(ContinuousFunction(domain),  [  self.survey.getStationLocation(ip) for ip in self.survey.getInjectionStations() ])

        self.field_resolution=field_resolution
        self.field_origin=field_origin
        self.sigma_background=sigma_background
        self.gamma_background=gamma_background
        self.padding_tags=padding_tags

        self.injections= [ i for i in self.survey.injectionIterator()]
        self.injectionMap=[ k for k in range(len(self.injections)) ]
        
        self.setUpDataMaps()
        self.setPrimaryPotential()
    

    
    def getAllInjections(self):
        return self.injectionMap
    
    def getInjection(self, k):
        return self.injections[self.injectionMap[k]]
    
    def getNumberOfInjections(self):
        return len(self.injections)
    
    def setUpDataMaps(self):
        """
        This sets up the mapping of the DC self.dataDCMaps[self.numSrc] and IP self.dataIPMaps[self.numSrc] predictions to an array d[self.numDataMax, self.numSrc]
        """
        self.numSrc=self.getNumberOfInjections()
        self.dataDCMaps={}
        self.dataIPMaps={}
        self.numData={}
        for k, i in enumerate(self.getAllInjections()):
                self.dataDCMaps[i] = { s: j for j,s in enumerate(self.survey.getObservations(self.getInjection(i)))}
                self.dataIPMaps[i] = { s: j+len(self.dataDCMaps[i]) for j,s in enumerate(self.survey.getObservations(self.getInjection(i)))}
                self.numData[i]=len(self.dataDCMaps[i])+len(self.dataIPMaps[i])
        self.numDataMax=max(self.numData.values())

        self.use=np.zeros((self.numDataMax, self.numSrc), dtype=bool)
        for k, i in enumerate(self.getAllInjections()):
                for j in self.dataDCMaps[i].values():
                    self.use[j,k]=True 
                for i in self.dataIPMaps[i].values():
                    self.use[j,k]=True
                    
    def makeDataSet(self, sources):
        """
        
        """
        responses=np.zeros((self.numDataMax, len(sources)), dtype=float)
        if self.survey.hasDipoleInjections():
            for k, ip in enumerate(sources):
                for s,i in self.dataDCMaps[ip].items():
                    responses[i,k]=self.survey.getDataRecord(self.getInjection(ip)+ s, datatype='R')
                for s,i in self.dataIPMaps[ip].items():
                    d=self.survey.getDataRecord( self.getInjection(ip) + s, datatype='R')
                    e=self.survey.getDataRecord(self.getInjection(ip) + s, datatype='ETA')
                    responses[i,k]=e/(1-e)*d
        else:
            for k, ip in enumerate(sources):
                for s,i in self.dataDCMaps[ip].items():
                    responses[i,k]=self.survey.getDataRecord( (self.getInjection(ip),) + s, datatype='R')
                for s,i in self.dataIPMaps[ip].items():
                    d=self.survey.getDataRecord( (self.getInjection(ip),) + s , datatype='R')
                    e=self.survey.getDataRecord( (self.getInjection(ip),) + s, datatype='ETA')
                    responses[i,k]=e/(1-e)*d            
        return responses
        
    def setPrimaryPotential(self):
        """
        this sets the primary potential assuming sigma=1 and I=1
        """
        self.primary_potential={}
        self.primary_potential_at_stations = {}
        self.pde.setValue(A=kronecker(3), X=Data())       
        for i, ip in enumerate(self.survey.getListOfInjectionStations()):
            s=Scalar(0.,DiracDeltaFunctions(self.domain))
            if self.stationsFMT is None:
                s.setTaggedValue(ip,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%ip,1.)
            self.pde.setValue(y_dirac=s)
            self.primary_potential[ip]=self.pde.getSolution()
            self.primary_potential_at_stations[ip]=np.array(self.observation_locator(self.primary_potential[ip]))
            print("Primary potential for %s: %s"%(ip,str(self.primary_potential[ip])))

    def runSurvey(self, sources, sigma_field, gamma_field):
        # sources point into 
        # array to return data: 
        responses=np.zeros((self.numDataMax, len(sources)), dtype=float)
        
        # extend the fields to the domain and grep values at source locations: 
        sigma, sigma_p=mapToDomain(self.domain, sigma_field, self.field_resolution, origin=self.field_origin, data0=self.sigma_background, tags0=self.padding_tags, locators=self.source_locator )
        gamma, gamma_p=mapToDomain(self.domain, gamma_field, self.field_resolution, origin=self.field_origin, data0=self.gamma_background, tags0=self.padding_tags, locators=self.source_locator )
    
        self.pde.setValue(A=sigma*kronecker(3), y_dirac=Data())
        secondary_potential_at_stations={}
        u_at_stations={}
        # DC .... 
        for k, j in enumerate(sources):
            if self.survey.hasDipoleInjections():
                ips=self.getInjection(j)
                for ip in ips:
                    if not ip in secondary_potential_at_stations:
                        idx=self.survey.getInjectionStationIndex(ip)
                        sigma0=sigma_p[idx]
                        print("DC injection %s at %s, sigma_p=%e"%(ip, idx, sigma0))

                        self.pde.setValue(X=(1-sigma/sigma0)*grad(self.primary_potential[ip])) 
                        u_s=self.pde.getSolution()
                        secondary_potential_at_stations[ip]=np.array(self.observation_locator(u_s))

                        u_at_stations[ip]=secondary_potential_at_stations[ip]+self.primary_potential_at_stations[ip]/sigma0            
                for s,i in self.dataDCMaps[j].items():
                    Midx, Nidx=self.survey.getObservationElectrodeIndex(s[0]), self.survey.getObservationElectrodeIndex(s[1]) 
                    responses[i,k]=u_at_stations[ips[0]][Midx]-u_at_stations[ips[0]][Nidx]- u_at_stations[ips[1]][Midx]+u_at_stations[ips[1]][Nidx]                       
            else:
                ip=self.getInjection(j)
                idx=self.survey.getInjectionStationIndex(ip)
                sigma0=sigma_p[idx]
                print("DC injection %s at %s, sigma_p=%e"%(ip, idx, sigma0))

                self.pde.setValue(X=(1-sigma/sigma0)*grad(self.primary_potential[ip])) 
                u_s=self.pde.getSolution()
                secondary_potential_at_stations[ip]=np.array(self.observation_locator(u_s))

                u_at_stations=secondary_potential_at_stations[ip]+self.primary_potential_at_stations[ip]/sigma0            
                for s,i in self.dataDCMaps[j].items():
                    Midx, Nidx=self.survey.getObservationElectrodeIndex(s[0]), self.survey.getObservationElectrodeIndex(s[1]) 
                    responses[i,k]=u_at_stations[Midx]-u_at_stations[Nidx]
        
        #.. IP
        sigma2=sigma/(1+gamma)
        du_at_stations={}
        u_s={}
        self.pde.setValue(A=sigma2*kronecker(3), y_dirac=Data())
        for k, j in enumerate(sources):
            
            if self.survey.hasDipoleInjections():
                ips=self.getInjection(j)
                for ip in ips:
                    if not ip in u_s:
                        idx=self.survey.getInjectionStationIndex(ip)
                        sigma20=sigma_p[idx]/(1+gamma_p[idx])
                        sigma0=sigma_p[idx]
                        print("IP injection %s at %s, sigma2_p, gamma_p = %e, %e"%(ip, idx, sigma20, gamma_p[idx]))
                        self.pde.setValue(X=(1-sigma2/sigma20)*grad(self.primary_potential[ip]))   
                        
                        u_s[ip]=self.pde.getSolution()
                        du_at_stations[ip]=np.array(self.observation_locator(u_s[ip]))-secondary_potential_at_stations[ip]+self.primary_potential_at_stations[ip]*(gamma_p[idx]/sigma0)
                for s,i in self.dataIPMaps[j].items():
                    Midx, Nidx=self.survey.getObservationElectrodeIndex(s[0]), self.survey.getObservationElectrodeIndex(s[1]) 
                    responses[i,k]=du_at_stations[ips[0]][Midx]-du_at_stations[ips[0]][Nidx]-du_at_stations[ips[1]][Midx]+du_at_stations[ips[1]][Nidx]
            else:
                ip=self.getInjection(j)
                idx=self.survey.getInjectionStationIndex(ip)
                sigma20=sigma_p[idx]/(1+gamma_p[idx])
                sigma0=sigma_p[idx]
                print("IP injection %s at %s, sigma2_p, gamma_p = %e, %e"%(ip, idx, sigma20, gamma_p[idx]))
                self.pde.setValue(X=(1-sigma2/sigma20)*grad(self.primary_potential[ip])) 
            
                u_s=self.pde.getSolution()
                du_at_stations=np.array(self.observation_locator(u_s))-secondary_potential_at_stations[ip]+self.primary_potential_at_stations[ip]*(gamma_p[idx]/sigma0)
                for s,i in self.dataIPMaps[j].items():
                    Midx, Nidx=self.survey.getObservationElectrodeIndex(s[0]), self.survey.getObservationElectrodeIndex(s[1]) 
                    responses[i,k]=du_at_stations[Midx]-du_at_stations[Nidx]
        
        self.sigma=sigma
        self.gamma=gamma
        
        return responses # [self.numDataMax, len(sources)]
        
