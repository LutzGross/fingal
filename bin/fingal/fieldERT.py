"""
tools for ERT inversion including ERT cost functions for inversion

by l.gross@uq.edu.au, 2021
"""

from esys.escript import *
import numpy as np
from esys.downunder import MeteredCostFunction
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag, getInfLocator
import logging
from esys.weipa import saveVTK, saveSilo

LOGCLIP=15


class FieldIntensityERTInversion(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., sigma0=.001, region_fixed=Data(), stationsFMT="e%s", 
                 weightLogDefect=0., adjustStationLocationsToElementCenter=True, logclip=15):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(FieldIntensityERTInversion, self).__init__()
        self.datatol=1e-30
        self.sigma0=sigma0
        self.stationsFMT=stationsFMT
        self.useLogDefect=useLogDefect

        if self.useLogDefect:
            lslogger.info("Misfit is using logarithm.") 
        else:
            lslogger.info("Misfit is using norm relative difference.")
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=Function(domain).getX()
        
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,X.getFunctionSpace())
            self.misfit.data[AB]=Scalar(0.,X.getFunctionSpace()) 
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, X.getFunctionSpace()))   
            for A,B in self.data.getInjections(M):
                D_ABS=abs(self.data.getFieldIntensityData((A,B,M)))+self.datatol
                self.misfit.w[(A,B)].copyWithMask(Scalar(1,self.misfit.w[AB].getFunctionSpace()), mask)   # 1 where data are measured @M
                self.misfit.data[(A,B)].copyWithMask(Scalar(D_ABS,self.misfit.w[AB].getFunctionSpace()), mask) # data inserted  @ M 

#                self.misfit.w[(A,B)]+=mask*(1-self.misfit.w[(A,B)])   # 1 where data are measured @M
#                self.misfit.data[(A,B)]+=mask*(D_ABS-self.misfit.data[(A,B)]) # data inserted  @ M 
        lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            s=integrate(self.misfit.w[AB])
            #self.misfit.data[AB]+=(1-wherePositive(self.misfit.w[AB])) # one inserted to avoid division by zero in misfit            
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # data inserted  @ M 

            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        lslogger.info("maximal station adjustment is %e"%adjustmax)
        lslogger.info("building primary electric fields (this will take some time)")
        self.phi_p=self.getPrimaryElectricPotentials(sigma0)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain, poisson=(abs(w1)>0) ) 
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain)        
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
            if not self.alpha0 > 0: 
                self.Spde.setValue(q=region_fixed)
        else:
            self.Spde=None
 
    def getPrimaryElectricPotentials(self, sigma):
        """
        return the primary electric for the injections (A,B)
        """
        primary_potential={}
        self.pde.setValue(A=sigma*kronecker(3), y_dirac=Data(), Y=Data(), X=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            primary_potential[A]=self.pde.getSolution()
            lslogger.debug("primary potential for injection at %d -> %s"%(A,str(primary_potential[A])))
        return primary_potential

    def getSecondaryElectricPotentials(self, sigma, sigma0, primary_potentials):
        secondary_potentials={}
        print("getSecondaryElectricPotentials: sigma=",str(sigma))
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
        return secondary_potentials
    
    def getSigma(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m*self.alpha0
        else:
            p=m
        return self.sigma0*exp(p)
    
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)

    def _getArguments(self, m):
        
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)

        sigma=self.getSigma(ppi, isSmoothed=True)
        secondary_potentials=self.getSecondaryElectricPotentials(sigma, self.sigma0, self.phi_p)

        return secondary_potentials, sigma, ppi           
    
    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        ppi=args[2]

        mi=interpolate(m, Function(self.pde.getDomain()))
        
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit=None
        for A,B in self.data.injectionIterator():
                if misfit is None:
                    misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], self.misfit.w[(A,B)].getFunctionSpace())
                if self.useLogDefect:
                    EI_AB=length(E_AB)+self.datatol
                    misfit+=self.misfit.w[(A,B)]*log(EI_AB/self.misfit.data[(A,B)])**2
                else:
                    EI_AB=length(E_AB)
                    misfit+=self.misfit.w[(A,B)] * (1-(EI_AB/self.misfit.data[(A,B)]))**2
        A2=integrate(misfit)
        lslogger.info("sigma = %s"%(str(sigma)))
        lslogger.debug("p = %s"%(str(ppi)))
        lslogger.debug("m = %s"%(str(m)))
        lslogger.debug("L2, H1, misfit= %e, %e, %e"%(A0/2, A1/2, A2/2))
        return (A0+A1+A2)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*interpolate(m, X.getFunctionSpace())
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            EI_AB=length(E_AB)+self.datatol
            if self.useLogDefect:       
                self.pde.setValue(X=self.misfit.w[(A,B)]/(EI_AB**2)*log(EI_AB/self.misfit.data[(A,B)])*E_AB )
            else:
                raise NotImplementedError
            ustar=self.pde.getSolution()
            Y2-=inner(grad(ustar, E_AB.getFunctionSpace()),E_AB)

        if self.Spde:
            self.Spde.setValue(Y=Y2*sigma)
            Y+=self.Spde.getSolution()
        else:
            Y+=Y2*sigma/self.alpha0
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        self.Hpde.setValue(X=r[1], Y=r[0])
        #saveVTK("test", PPP=self.Hpde.getRightHandSide())
        p=self.Hpde.getSolution()
        #p*=1./Lsup(p)
        lslogger.debug("inverse Hessian called. search direction = %s",p)
        return p
    
class ChargeabilityLinearizedInversion(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0.,  gamma0=0.001, sigma=0.001, region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, useLogDefect=True):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(ChargeabilityLinearizedInversion, self).__init__()
        self.datatol=1e-30
        self.useLogDefect=useLogDefect
        self.sigma=sigma
        self.gamma0=gamma0
        self.useDifferenceOfFields=False
        self.stationsFMT=stationsFMT
        self.misfitFunctionSpace=Function(domain)
        if getMPIRankWorld() == 0:
            lslogger.debug("Costfunction = (E_relaxed-E).E/|E|^2")  
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))
        self.pde.setValue(A=self.sigma*kronecker(3))
         
        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        if getMPIRankWorld() == 0:
             lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=self.misfitFunctionSpace.getX()
        #X=ReducedFunction(domain).getX()
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,self.misfitFunctionSpace)
            self.misfit.data[AB]=Scalar(0.,self.misfitFunctionSpace)
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
                
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, self.misfitFunctionSpace))   
            for A,B in self.data.getInjections(M):
                D_ABM=self.data.getChargeabilityData((A,B,M))
                E_ABM=self.data.getFieldIntensityData((A,B,M))
                if abs(D_ABM) > 0.:
                    self.misfit.w[(A,B)].copyWithMask(Scalar(1./D_ABM**2,self.misfitFunctionSpace), mask)   # 1 where data are measured @M
                    self.misfit.data[(A,B)].copyWithMask(Scalar(D_ABM,self.misfitFunctionSpace), mask) # data inserted  @ M 
                    
        if adjustStationLocationsToElementCenter and getMPIRankWorld() == 0:
            lslogger.info("maximal station adjustment is %e"%adjustmax)
        if getMPIRankWorld() == 0:
            lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # insert 1's to avoid division by zero
            s=integrate(self.misfit.w[AB]*self.misfit.data[AB]**2)
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        self.phi_p=self.getElectricPotentials()
        #self.secondary_potentials=self.getSecondaryElectricPotentials(self.sigma, self.sigma0, self.phi_p)
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain)        
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        
        self.Spde=None
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain) 
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
            if not self.alpha0>0:
                    self.Spde.setValue(q=region_fixed)
 
    def getElectricPotentials(self):
        """
        return the primary electric for the injections (A,B)
        """
        potential={}
        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            potential[A]=self.pde.getSolution()
        if getMPIRankWorld() == 0:
            lslogger.debug("primary potential for %s injection calculated."%len(potential))
        return potential

    def getSecondaryElectricPotentials(self, gamma, primary_potentials):
        secondary_potentials={}
        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=gamma*self.sigma*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
        if getMPIRankWorld() == 0:
            lslogger.debug("%s secondary potential calculated"%(len(secondary_potentials)))
        return secondary_potentials
            
    def getCargeability(self, m, isSmoothed=False):
        """
        return Cargeability (eta) for a given property function m
        """
        gamma=self.getGamma(m, isSmoothed)
        return gamma/(1.+gamma)
    
    def getGamma(self, m, isSmoothed=False):
        """
        
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m/self.alpha0
        else:
            p=m
        gamma=self.gamma0*exp(clip(p, minval=-LOGCLIP, maxval=LOGCLIP))
        return gamma
        
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)


    def _getArguments(self, m):
        if self.Spde :
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
        else:
            p=m/self.alpha0
                
        gamma=self.getGamma(p, isSmoothed=True)
        gammai=interpolate(gamma,  Function(self.pde.getDomain()))
        secondary_potentials=self.getSecondaryElectricPotentials(gammai, self.phi_p)
        return gammai, secondary_potentials, 

    #==snip ===============================================
    
    
    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        gammai=args[0]
        secondary_potentials=args[1]

        mi=interpolate(m,  Function(self.pde.getDomain()))
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        
        misfit=Scalar(0.,self.misfitFunctionSpace )
        for A, B in self.data.injectionIterator():
            E =-grad(self.phi_p[A]-self.phi_p[B], misfit.getFunctionSpace())
            DE=-grad(secondary_potentials[A]-secondary_potentials[B], misfit.getFunctionSpace())
            L_E2=length(E)**2
            #misfit+=self.misfit.w[(A,B)]*safeDiv(L_E2*self.misfit.data[(A,B)] - inner(DE, E), L_E2*self.misfit.data[(A,B)])**2
            misfit+=self.misfit.w[(A,B)]*((L_E2*self.misfit.data[(A,B)] - inner(DE, E))/(L_E2+self.datatol**2))**2
        A2=integrate(misfit)

        if lslogger.isEnabledFor(logging.DEBUG):
            strgamma=str(gammai)
            strm=str(mi)
            if getMPIRankWorld() == 0:
                lslogger.debug("gamma = %s"%strgamma)
                lslogger.debug("m = %s"%strm)
                lslogger.debug("L2, H1, misfit = %e, %e, %e"%(A0/2, A1/2, A2/2))

        return (A0+A1+A2)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        gammai=args[0]
        secondary_potentials=args[1]
        mi=interpolate(m,  Function(self.pde.getDomain()))

        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*mi

        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        Y2=Scalar(0.,self.misfitFunctionSpace )
        for A, B in self.data.injectionIterator():
            E =-grad(self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            DE=-grad(secondary_potentials[A]-secondary_potentials[B], Y2.getFunctionSpace())
            L_E2=length(E)**2
            #self.pde.setValue(X=self.misfit.w[(A,B)]*safeDiv(self.misfit.data[(A,B)]*L_E2 - inner(DE, E), (self.misfit.data[(A,B)]*L_E2)**2)* E)
            self.pde.setValue(X=self.misfit.w[(A,B)]*(self.misfit.data[(A,B)]*L_E2 - inner(DE, E))/(L_E2+self.datatol**2)**2* E)
            ustar=self.pde.getSolution()
            Y2+=-inner(grad(ustar, E.getFunctionSpace()),E)
        if self.Spde:
            self.Spde.setValue(Y=Y2*gammai*self.sigma)
            Y+=self.Spde.getSolution()
        else:                
            Y+=Y2*gammai*self.sigma/self.alpha0

        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        mi=interpolate(m,  Function(self.pde.getDomain()))
        if getMPIRankWorld() == 0:
            lslogger.debug("inverse Hessian called.")
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        return p

class ChargeabilityLinearizedInversion2(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    POSEXP=2
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., sigma=0.001, beta=10., region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, useLogDefect=True):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(ChargeabilityLinearizedInversion, self).__init__()
        self.NONLINEAR2=False
        self.datatol=1e-30
        self.useLogDefect=useLogDefect
        self.sigma=sigma
        self.useDifferenceOfFields=False
        self.stationsFMT=stationsFMT
        if self.useLogDefect:
            raise ValueError("Log Defect not supported yet.")
            if self.useDifferenceOfFields:
                lslogger.info("Costfunction = log((E_relaxed-E).E_relaxed/|E_relaxed|^2)") 
            else:
                lslogger.info("Costfunction = log((E_relaxed-E).E/(E_relaxed.E))")
        else:
            if self.NONLINEAR2:
                lslogger.info("Costfunction = (E_relaxed-E).E_relaxed/|E_relaxed|^2") 

            else:
                lslogger.info("Costfunction = (E_relaxed-E).E/|E|^2") 
            
            
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))
        self.pde.setValue(A=self.sigma*kronecker(3))
         
        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=Function(domain).getX()
        #X=ReducedFunction(domain).getX()
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,X.getFunctionSpace())
            self.misfit.data[AB]=Scalar(0.,X.getFunctionSpace())
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
                
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, X.getFunctionSpace()))   
            for A,B in self.data.getInjections(M):
                D_ABM=self.data.getChargeabilityData((A,B,M))
                E_ABM=self.data.getFieldIntensityData((A,B,M))
                if abs(D_ABM) > 0.:
                    self.misfit.w[(A,B)].copyWithMask(Scalar(1,self.misfit.w[AB].getFunctionSpace()), mask)   # 1 where data are measured @M
                    self.misfit.data[(A,B)].copyWithMask(Scalar(D_ABM,self.misfit.w[AB].getFunctionSpace()), mask) # data inserted  @ M 
        if adjustStationLocationsToElementCenter:
            lslogger.info("maximal station adjustment is %e"%adjustmax)
        lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            s=integrate(self.misfit.w[AB]*self.misfit.data[AB]**2)
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # insert 1's to avoid division by zero
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:

        lslogger.info("building primary electric fields (this will take some time)")
        self.phi_p=self.getPrimaryElectricPotentials()
        #self.secondary_potentials=self.getSecondaryElectricPotentials(self.sigma, self.sigma0, self.phi_p)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        self.beta=beta#/vol(domain)
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain)        
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        
        self.Spde=None
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain) 
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
            if not self.alpha0>0:
                    self.Spde.setValue(q=region_fixed)
 
    def getPrimaryElectricPotentials(self):
        """
        return the primary electric for the injections (A,B)
        """
        primary_potential={}
        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            primary_potential[A]=self.pde.getSolution()
            lslogger.debug("primary potential for injection at %d -> %s"%(A,str(primary_potential[A])))
        return primary_potential

    def getSecondaryElectricPotentials(self, gamma, primary_potentials):
        secondary_potentials={}
        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=gamma*self.sigma*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
        lslogger.debug("%s secondary potential for injections calculated"%(len(secondary_potentials)))
        return secondary_potentials
            
    def getSigmaIP(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                gamma=self.Spde.getSolution()
            else:
                gamma=m/self.alpha0
        else:
            gamma=m

        return self.sigma*(1-gamma)
    
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)

    def _getArguments(self, m):
        
        if self.Spde:
            self.Spde.setValue(Y=m)
            gamma=self.Spde.getSolution()
        else:
            gamma=m/self.alpha0
        gammai=interpolate(gamma, Function(self.pde.getDomain()))
            
        secondary_potentials=self.getSecondaryElectricPotentials(gammai, self.phi_p)
        return secondary_potentials, gammai



    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        gammai=args[1]
        mi=interpolate(m,  Function(self.pde.getDomain()))

        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit=None
        for A, B in self.data.injectionIterator():
            if misfit is None:
                misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                

            E =-grad(self.phi_p[A]-self.phi_p[B], misfit.getFunctionSpace())
            DE=-grad(secondary_potentials[A]-secondary_potentials[B])
            
            if self.NONLINEAR2:
                L_EHAT2=length(E+DE)**2
                ETA=safeDiv(inner(DE, E+DE), L_EHAT2)
            else:
                L_E2=length(E)**2
                ETA=safeDiv(inner(DE, E), L_E2)
            #MIS2=(1- ETA / self.misfit.data[(A,B)])**2
            MIS2=(self.misfit.data[(A,B)]- ETA)**2
            misfit+=self.misfit.w[(A,B)]*MIS2
        A2=integrate(misfit)
        t=gammai*(gammai-1)
        A3=integrate((t*wherePositive(t))**self.POSEXP)*self.beta/self.POSEXP
        if getMPIRankWorld() == 0:
            lslogger.debug("m = %s"%(str(mi)))
            lslogger.debug("gamma = %s"%(str(gammai)))
            lslogger.debug("L2, H1, misfit, positive= %e, %e, %e, %e"%(A0/2, A1/2, A2/2, A3/2))
        return (A0+A1+A2+A3)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        gammai=args[1]
        mi=interpolate(m,  Function(self.pde.getDomain()))

        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*mi #+self.beta**self.POSEXP * (mi*wherePositive(mi))**(self.POSEXP-1)
        self.pde.setValue(X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E =-grad(self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            DE=-grad(secondary_potentials[A]-secondary_potentials[B], Y2.getFunctionSpace())
            
            if self.NONLINEAR2:
                EHAT=E+DE
                L_EHAT2=length(EHAT)**2
                ETA=safeDiv(inner(DE, EHAT), L_EHAT2)
                MIS=self.misfit.data[(A,B)]- ETA 
                self.pde.setValue(X=self.misfit.w[(A,B)]*MIS*(2 * inner(EHAT,E)/L_EHAT2**2*EHAT-E/L_EHAT2))
            else:
                L_E2=length(E)**2
                ETA=safeDiv(inner(DE, E), L_E2)
                MIS=self.misfit.data[(A,B)]- ETA  
                self.pde.setValue(X=self.misfit.w[(A,B)]*MIS/L_E2*E)
            ustar=self.pde.getSolution()
            Y2+=-inner(grad(ustar, E.getFunctionSpace()),E)
            
        t=gammai*(gammai-1)
        POS=wherePositive(t)*(t*wherePositive(t))**(self.POSEXP-1)*self.beta*(gammai-0.5)
        if self.Spde:
            if Y2.getFunctionSpace() == ReducedFunction(self.Spde.getDomain()):
                self.Spde.setValue(Y=POS, Y_reduced=Y2*self.sigma)
            else:
                self.Spde.setValue(Y=POS+Y2*self.sigma)
            Y+=self.Spde.getSolution()
        else:                
            Y+=(POS+Y2*self.sigma)/self.alpha0
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        mi=interpolate(m,  Function(self.pde.getDomain()))
        lslogger.debug("inverse Hessian called.")
        #if self.beta >0 and self.POSEXP>0:
        #    t=mi*(mi-1)
        #    tp=wherePositive(t)
        #    if self.POSEXP ==1:
        #        D=tp 
        #    elif self.POSEXP ==2:
        #        D=tp*(t +t**2*(mi-0.5)**2) 
        #    else:
        #        D=tp*((t*tp)**(self.POSEXP-1) +(t*tp)**2*(self.POSEXP-1)*(mi-0.5)**2*(t*tp)**(self.POSEXP-2)) 

        #else: 
        #    D=0
        #print( "D_beta =%s"%(str(D)) )
        #self.Hpde.setValue(D=self.w0+D*self.beta, X=r[1], Y=r[0])
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        return p

class ChargeabilityInversion(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    CASE= 3
    POSEXP=3
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., sigma=0.001, sigma0=0.001, beta=10., region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, useLogDefect=True, useDifferenceOfFields=False, sigmatest=None):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(ChargeabilityInversion, self).__init__()
        self.datatol=1e-30
        self.useLogDefect=useLogDefect
        self.sigma0=sigma0
        self.sigma=sigma
        self.stationsFMT=stationsFMT
        self.useDifferenceOfFields=useDifferenceOfFields
        if sigmatest is None:
            self.sigmatest=sigma
        else:
            self.sigmatest=sigmatest
        if self.useLogDefect:
            if self.useDifferenceOfFields:
                lslogger.info("Costfunction = log((E_relaxed-E).E_relaxed/|E_relaxed|^2)") 
            else:
                lslogger.info("Costfunction = log((E_relaxed-E).E/(E_relaxed.E))")
        else:
            if self.useDifferenceOfFields:
                lslogger.info("Costfunction = (E_relaxed-E).E_relaxed/|E_relaxed|^2") 
            else:
                lslogger.info("Costfunction = (E_relaxed-E).E/(E_relaxed.E)")
        
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()
        if adjustStationLocationsToElementCenter:
            lslogger.info("maximal station adjustment is %e"%adjustmax)
        ########################################################
        lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=Function(domain).getX()
        #X=ReducedFunction(domain).getX()
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,X.getFunctionSpace())
            self.misfit.data[AB]=Scalar(0.,X.getFunctionSpace())
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
                
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, X.getFunctionSpace()))   
            for A,B in self.data.getInjections(M):
                D_ABM=self.data.getChargeabilityData((A,B,M))
                E_ABM=self.data.getFieldIntensityData((A,B,M))
                if abs(D_ABM) > 0.:
                    if self.CASE in [3, 5]:
                        #self.misfit.w[(A,B)]+=mask*(1.-self.misfit.w[(A,B)])   # =1 at sensor locations M,  0 elsewhere
                        #self.misfit.data[(A,B)]+=mask*(D_ABM-self.misfit.data[(A,B)]) # = data at sensor locations M
                        
                        self.misfit.w[(A,B)].copyWithMask(Scalar(1,self.misfit.w[AB].getFunctionSpace()), mask)   # 1 where data are measured @M
                        self.misfit.data[(A,B)].copyWithMask(Scalar(D_ABM,self.misfit.w[AB].getFunctionSpace()), mask) # data inserted  @ M 

                    else: # obsolete
                        if self.CASE in [8, 9, 10, 11,12]:
                            self.misfit.w[(A,B)]+=mask*(1.-self.misfit.w[(A,B)])
                            self.misfit.data[(A,B)]+=mask*(D_ABM-self.misfit.data[(A,B)])
                        else:
                            self.misfit.w[(A,B)]+=mask*(1/D_ABM**2-self.misfit.w[(A,B)])
                            self.misfit.data[(A,B)]+=mask*(D_ABM-self.misfit.data[(A,B)])

        lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            if self.CASE in [3, 5]:
                s=integrate(self.misfit.w[AB])
                #self.misfit.data[AB]+=(1-wherePositive(self.misfit.w[AB])) # insert 1's to avoid division by zero
                self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # insert 1's to avoid division by zero
            else:
                s=integrate(self.misfit.w[AB]*self.misfit.data[AB]**2)
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        lslogger.info("building primary electric fields (this will take some time)")
        self.phi_p=self.getPrimaryElectricPotentials(self.sigma0)
        self.secondary_potentials=self.getSecondaryElectricPotentials(self.sigma, self.sigma0, self.phi_p)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        self.beta=beta#/vol(domain)
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain)        
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain)        
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0) #, q=region_fixed)
        else:
            self.Spde=None
 
    def getPrimaryElectricPotentials(self, sigma ):
        """
        return the primary electric for the injections (A,B)
        """
        primary_potential={}
        self.pde.setValue(A=sigma*kronecker(3), y_dirac=Data(), Y=Data(), X=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            primary_potential[A]=self.pde.getSolution()
            lslogger.debug("primary potential for injection at %d -> %s"%(A,str(primary_potential[A])))
        return primary_potential

    def getSecondaryElectricPotentials(self, sigma, sigma0, primary_potentials):
        secondary_potentials={}
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
            #lslogger.debug("secondary potential for injection at %d -> %s"%(A,str(secondary_potentials[A])))
        return secondary_potentials
            
    def getSigmaIP(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m*self.alpha0
        else:
            p=m

        return self.sigmatest*exp(p)
    
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)

    def _getArguments(self, m):
        
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)
        else:
            ppi=clip(interpolate(m*self.alpha0, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)

        relaxedsigma=self.getSigmaIP(ppi, isSmoothed=True)
        secondary_potentials_relaxed=self.getSecondaryElectricPotentials(relaxedsigma, self.sigma0, self.phi_p)
        return secondary_potentials_relaxed, relaxedsigma, ppi



    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials_relaxed=args[0]
        relaxedsigma=args[1]
        ppi=args[2]
        mi=interpolate(m, ppi.getFunctionSpace())
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit=None
        for A, B in self.data.injectionIterator():
            if misfit is None:
                misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                
            E_AB =-grad(self.secondary_potentials[A]-self.secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], misfit.getFunctionSpace())
            DE_AB=-grad(secondary_potentials_relaxed[A]-secondary_potentials_relaxed[B] -self.secondary_potentials[A]+self.secondary_potentials[B], misfit.getFunctionSpace())

            DxE=inner(DE_AB, E_AB)
            E_AB_relaxed=E_AB+DE_AB
            L_E_AB_relaxed=length(E_AB_relaxed)
            #R_AB=safeDiv(D, self.misfit.data[(A,B)] )
            #E_AB_relaxed=E_AB+DE_AB

            #R_AB=safeDiv(inner(DE_AB, E_AB), inner(E_AB_relaxed, E_AB)*self.misfit.data[(A,B)] )            
            # E_AB_relaxed=-grad(secondary_potentials_relaxed[A]-secondary_potentials_relaxed[B]+self.phi_p[A]-self.phi_p[B], misfit.getFunctionSpace())
            
            
            #f=self.misfit.w[(A,B)]/Lsup(self.misfit.w[(A,B)])
            #print f*eta_AB/self.misfit.data[(A,B)], f*eta_AB, f*self.misfit.data[(A,B)]
            
            if self.CASE in [3,5]:
                if self.useDifferenceOfFields: # CASE 5
                    ETA=safeDiv(inner(DE_AB, E_AB_relaxed), length(L_E_AB_relaxed)**2)
                else: # CASE==3 
                    ETA=safeDiv(inner(DE_AB, E_AB),inner(E_AB_relaxed, E_AB))
                if self.useLogDefect:
                    MIS=log( self.datatol+(abs(ETA/self.misfit.data[(A,B)])) )**2
                else:
                    MIS=(1- ETA / self.misfit.data[(A,B)])**2


            else:    # obsolete bits:
                if self.CASE==1:
                    ETA=safeDiv(length(DE_AB),length(E_AB+DE_AB))
                elif self.CASE==2:
                    ETA=safeDiv(DxE+length(DE_AB)**2,DxE+length(E_AB)**2)
                elif self.CASE==4:
                    ETA=DxE+length(DE_AB)**2
                elif self.CASE==6:
                    ETA=length(DE_AB)
                elif self.CASE==7:
                    ETA=length(E_AB+DE_AB)
                elif self.CASE==8:
                    ETA=1-length(E_AB)/(length(E_AB+DE_AB)+self.datatol)
                elif self.CASE==9:
                    ETA=DxE/(length(E_AB)+self.datatol)
                elif self.CASE==10:
                    ETA=DxE/(length(E_AB)+self.datatol)**2
                elif self.CASE==11:
                    ETA=log(length(E_AB+DE_AB)/(length(E_AB)+self.datatol)**2+self.datatol)
                elif self.CASE==12:
                    ETA=log(length(DE_AB)+self.datatol)
                elif self.CASE==13:
                    ETA=log(length(E_AB+DE_AB)+self.datatol)
                else:
                    raise NotImplementedError
                MIS=(ETA-self.misfit.data[(A,B)])**2
            misfit+=self.misfit.w[(A,B)]*MIS
                
                
                #f=wherePositive(self.misfit.w[(A,B)])
                #print("DIff =",str(f*length(DE_AB)), str(f*length(E_AB)), str(f*length(E_AB+DE_AB)))
                #print("D", str(f*length(E_AB+DE_AB)/(length(E_AB)+self.datatol)**2), str(f*MIS), str(f*self.misfit.data[(A,B)]))
                #saveVTK("x", d=f*length(E_AB+DE_AB)/(length(E_AB)+self.datatol)**2, E=f*length(E_AB), E2=f*length(E_AB+DE_AB))
                #1/0
                #misfit+=self.misfit.w[(A,B)] * (1-(eta_AB/self.misfit.data[(A,B)]))**2
                #print("D =", str(f*(D-self.misfit.data[(A,B)])**2/(self.misfit.data[(A,B)]**2+self.datatol)))         
                #misfit+=self.misfit.w[(A,B)]*(D-self.misfit.data[(A,B)])**2/(self.misfit.data[(A,B)]**2+self.datatol)
                #print("eta(%s) ="%self.CASE, str(ETA*f), "misfit^2 =", str(MIS*f), "DE =",str(DE_AB*f))
                #print("total = ", integrate(self.misfit.w[(A,B)]*MIS), integrate(self.misfit.w[(A,B)]*self.misfit.data[(A,B)]**2))
        
        A2=integrate(misfit)
        A3=self.beta**self.POSEXP *integrate((ppi*wherePositive(ppi))**self.POSEXP)*2./self.POSEXP
        #A3=self.beta**self.POSEXP * integrate((mi*wherePositive(mi))**self.POSEXP)*2./self.POSEXP
        if getMPIRankWorld() == 0:
            lslogger.info("relaxedsigma = %s"%(str(relaxedsigma)))
            lslogger.debug("p = %s"%(str(ppi)))
            lslogger.debug("m = %s"%(str(mi)))
            lslogger.debug("L2, H1, misfit, positive= %e, %e, %e, %e"%(A0/2, A1/2, A2/2, A3/2))
        return (A0+A1+A2+A3)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials_relaxed=args[0]
        relaxedsigma=args[1]
        ppi=args[2]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        mi=interpolate(m, X.getFunctionSpace())
        Y=self.w0*mi #+self.beta**self.POSEXP * (mi*wherePositive(mi))**(self.POSEXP-1)
        self.pde.setValue(A=relaxedsigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB =-grad(self.secondary_potentials[A]-self.secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            DE_AB=-grad(secondary_potentials_relaxed[A]-secondary_potentials_relaxed[B] -self.secondary_potentials[A]+self.secondary_potentials[B], Y2.getFunctionSpace())
            #print "DE_AB", secondary_potentials_relaxed[A]-self.secondary_potentials[A], self.secondary_potentials[B]-secondary_potentials_relaxed[B]
            

            DxE=inner(DE_AB, E_AB)
            E_AB_relaxed=E_AB+DE_AB
            L_E_AB_relaxed=length(E_AB_relaxed)
            #DD=inner(E_AB_relaxed, E_AB)
            
            #D=inner(DE_AB, E_AB)
            #eta_AB=safeDiv(D, DD)
            #L_AB=length(E_AB)
            #    else: # CASE==3 
            
            if self.CASE in [3,5]:
                #if self.useDifferenceOfFields: # CASE 5
                #    #ETA=DxE
                #    ETA=safeDiv(inner(DE_AB, E_AB_relaxed), length(L_E_AB_relaxed)**2)
                #    MIS=(1-ETA/self.misfit.data[(A,B)])
                #    M=(safeDiv(1,length(L_E_AB_relaxed)**2)*E_AB-2*safeDiv(inner(E_AB_relaxed, E_AB),length(L_E_AB_relaxed)**4)*E_AB_relaxed)/self.misfit.data[(A,B)]
                #else: # CASE==3 
                #    ExE_relaxed=inner(E_AB_relaxed, E_AB)
                #    ETA=safeDiv(inner(DE_AB, E_AB),ExE_relaxed)
                #    MIS=(1.-ETA/self.misfit.data[(A,B)])
                #    F=safeDiv(length(E_AB)**2, ExE_relaxed**2)/self.misfit.data[(A,B)]
                #    M=-F*E_AB

                if self.useDifferenceOfFields: # CASE 5
                    ETA=safeDiv(inner(DE_AB, E_AB_relaxed), length(L_E_AB_relaxed)**2)
                    DETA=-(safeDiv(1,length(L_E_AB_relaxed)**2)*E_AB-2*safeDiv(inner(E_AB_relaxed, E_AB),length(L_E_AB_relaxed)**4)*E_AB_relaxed)
                else: # CASE==3 
                    ExE_relaxed=inner(E_AB_relaxed, E_AB)
                    ETA=safeDiv(inner(DE_AB, E_AB),ExE_relaxed)
                    DETA=safeDiv(length(E_AB)**2, ExE_relaxed**2)*E_AB
                if self.useLogDefect:
                    MIS=log(self.datatol+abs(ETA / self.misfit.data[(A,B)]) )
                    F=safeDiv(MIS,ETA)
                else:
                    MIS=(1- ETA / self.misfit.data[(A,B)])
                    F=-MIS/self.misfit.data[(A,B)]
                M=F*DETA

                #f=self.misfit.w[(A,B)]/Lsup(self.misfit.w[(A,B)])                
                #print "F=", F*f
                #print "DETA=", f*DETA
                #print "ETA=", f*DE_AB, f*E_AB_relaxed, f*length(L_E_AB_relaxed)**2
                #print "E_Ab_r=", E_AB_relaxed*f
                #print "MIS=", MIS*f
                #print "M=", M*f

            else:    # obsolete bits:
                if self.CASE==1:
                    ETA=safeDiv(length(DE_AB),length(E_AB_relaxed))
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=(DE_AB*length(E_AB_relaxed)**2-E_AB_relaxed*length(DE_AB)**2)/(length(E_AB_relaxed)**3*length(DE_AB)+self.datatol)
                elif self.CASE==4:
                    ETA=DxE+length(DxE)**2
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=2*E_AB_relaxed-E_AB
                elif self.CASE==6:
                    ETA=length(DE_AB)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=DE_AB/(ETA+self.datatol)
                elif self.CASE==7:
                    ETA=length(E_AB+DE_AB)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=(E_AB+DE_AB)/length(E_AB+DE_AB)
                elif self.CASE==8:
                    ETA=1-length(E_AB)/(length(E_AB+DE_AB)+self.datatol)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=E_AB_relaxed*length(E_AB)/(length(E_AB_relaxed)+self.datatol)**3
                elif self.CASE==9:
                    ETA=DxE/(length(E_AB)+self.datatol)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=E_AB/(length(E_AB)+self.datatol)
                elif self.CASE==10:
                    ETA=DxE/(length(E_AB)**2+self.datatol)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=E_AB/(length(E_AB)+self.datatol)**2
                elif self.CASE==11:
                    ETA=log(length(E_AB+DE_AB)/(length(E_AB)**2+self.datatol))
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=(E_AB+DE_AB)/(length(E_AB+DE_AB)+self.datatol)**2
                elif self.CASE==12:
                    ETA=ETA=log(length(DE_AB)+self.datatol)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=DE_AB/length(DE_AB)**2
                elif self.CASE==13:
                    ETA=log(length(E_AB_relaxed)+self.datatol)
                    MIS=(ETA-self.misfit.data[(A,B)])
                    M=E_AB_relaxed/(length(E_AB_relaxed)+self.datatol)**2
                else:
                    raise NotImplementedError
                M*=MIS
                #self.pde.setValue(X=self.misfit.w[(A,B)]/(self.misfit.data[(A,B)]*EI_AB_relaxed**2)*(1./self.misfit.data[(A,B)]-1./(eta_AB+self.datatol))*R_AB)
                #print "etaAB", eta_AB*f
                
                #print "E", E_AB_relaxed*f
                #print "R= ",R_AB/(length(DE_AB)*EI_AB_relaxed+self.datatol**2)*f
                
                #F=safeDiv((eta_AB-self.misfit.data[(A,B)])*L_AB**2,(self.misfit.data[(A,B)]*DD)**2) 
                #F=(eta_AB-self.misfit.data[(A,B)])*L_AB**2/((self.misfit.data[(A,B)]*DD)**2+self.datatol) 
                #F=(eta_AB-self.misfit.data[(A,B)])*L_AB**2/(DD**2+self.datatol)
                
                #F=(D-self.misfit.data[(A,B)])/(self.misfit.data[(A,B)]**2+self.datatol)
                #..F=(eta_AB-self.misfit.data[(A,B)])/(self.misfit.data[(A,B)]**2+self.datatol)*L_AB**2/(DD**2+self.datatol)
                #saveDataCSV("x.csv", F=F, N=(eta_AB-self.misfit.data[(A,B)])*L_AB**2, Z=(self.misfit.data[(A,B)]*DD)**2, DD=DD, e=self.misfit.data[(A,B)],eta_AB=eta_AB, L_AB=L_AB ,mask=f)
                #print "ETA=", (eta_AB-self.misfit.data[(A,B)])*f
                #print "DD=", DD*f
                #print "E_AB=", E_AB*f
                #print "L EAB=", f*L_AB**2
                #print "E_AB_relaxed=", E_AB_relaxed*f
                
                #print "z = ", f*L_AB**2/(DD**2+self.datatol)*E_AB[0]
                #print "z = ", f*L_AB**2/(DD**2+self.datatol)*E_AB[1]
                #print "z = ", f*L_AB**2/(DD**2+self.datatol)*E_AB[2]
                #print "D=", (eta_AB-self.misfit.data[(A,B)])*f
               
                
                #print("F0= ",str(F*E_AB[0]))
                #print("F1= ",str(F*E_AB[1]))
                ##print("F2= ",str(F*E_AB[2]))
                #1/0
            self.pde.setValue(X=self.misfit.w[(A,B)]*M)
            ustar=self.pde.getSolution()
            #print "ustar =", ustar*f
            #print "G*=", inner(grad(ustar, Y2.getFunctionSpace()),E_AB_relaxed)*relaxedsigma
            Y2-=inner(grad(ustar, Y2.getFunctionSpace()),E_AB_relaxed)
            #print "Y2=", Y2
#        saveSilo("x", mask=self.misfit.w[(A,B)]/Lsup(self.misfit.w[(A,B)])*log((eta_AB+self.datatol)/self.misfit.data[(A,B)]), E=self.misfit.w[(A,B)]/Lsup(self.misfit.w[(A,B)])*E_AB_relaxed)
        POS=self.beta**self.POSEXP * (ppi*wherePositive(ppi))**(self.POSEXP-1)
        if self.Spde:
            #print "Y start =",Y
            #print "Y2 =",Y2*relaxedsigma
            if Y2.getFunctionSpace() == ReducedFunction(self.Spde.getDomain()):
                self.Spde.setValue(Y=POS, Y_reduced=Y2*relaxedsigma)
            else:
                self.Spde.setValue(Y=POS+Y2*relaxedsigma)
            #print "Y2 bar =", self.Spde.getSolution()
            #
            #saveVTK("x", bar=self.Spde.getSolution(), Y=self.beta*ppi*wherePositive(ppi)+Y2*relaxedsigma, f=self.misfit.w[(A,B)]/Lsup(self.misfit.w[(A,B)]))
            Y+=self.Spde.getSolution()
        else:
                
            Y+=(POS+Y2*relaxedsigma)/self.alpha0
        #print "grad Y=",Y
        
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        #ppi=args[2]
        #mi=interpolate(m, r[0].getFunctionSpace())
        lslogger.debug("inverse Hessian called.")
        #self.Hpde.setValue(X=r[1], Y=r[0])
        #print "m = ", mi 
        #print "X,Y = ", r[1], r[0] 
        #self.Hpde.setValue(D=self.w0+self.beta**self.POSEXP*(self.POSEXP-1)*(mi*wherePositive(mi))**(self.POSEXP-2), X=r[1], Y=r[0])
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        #print "p = ", p 
        #saveVTK("x", p=p, Y=r[0], m=m)
        #1/0
        return p

class FieldIntensityInversionWithMinMean(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., penalty=1e-3, sigma0=.001, region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, useLogDefect=True):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(FieldIntensityInversionWithMinMean, self).__init__()
        self.datatol=1e-30
        self.sigma0=sigma0
        self.stationsFMT=stationsFMT
        self.useLogDefect=useLogDefect

        if self.useLogDefect:
            lslogger.info("Misfit is using logarithm.") 
        else:
            lslogger.info("Misfit is using norm relative difference.")
        # setup PDE:
        self.pde=setupERTPDE(domain)

        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=Function(domain).getX()
        
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,X.getFunctionSpace())
            self.misfit.data[AB]=Scalar(0.,X.getFunctionSpace()) 
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, X.getFunctionSpace()))   
            for A,B in self.data.getInjections(M):
                D_ABS=abs(self.data.getFieldIntensityData((A,B,M)))+self.datatol
                self.misfit.w[(A,B)].copyWithMask(Scalar(1,self.misfit.w[AB].getFunctionSpace()), mask)   # 1 where data are measured @M
                self.misfit.data[(A,B)].copyWithMask(Scalar(D_ABS,self.misfit.w[AB].getFunctionSpace()), mask) # data inserted  @ M 

#                self.misfit.w[(A,B)]+=mask*(1-self.misfit.w[(A,B)])   # 1 where data are measured @M
#                self.misfit.data[(A,B)]+=mask*(D_ABS-self.misfit.data[(A,B)]) # data inserted  @ M 
        lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            s=integrate(self.misfit.w[AB])
            #self.misfit.data[AB]+=(1-wherePositive(self.misfit.w[AB])) # one inserted to avoid division by zero in misfit            
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # data inserted  @ M 

            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        lslogger.info("maximal station adjustment is %e"%adjustmax)
        lslogger.info("building primary electric fields (this will take some time)")
        self.phi_p=self.getPrimaryElectricPotentials(sigma0)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain, poisson=(abs(w1)>0) ) 
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain)        
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
            if not self.alpha0 > 0: 
                self.Spde.setValue(q=region_fixed)
        else:
            self.Spde=None
        
        # this is bit preliminary:
        self.penalty=penalty
        self.u=Scalar(1., Solution(domain))
        self.u*=1./sqrt(integrate(self.w0*self.u*self.u))
        self.VOL=vol(domain)
        print("u =",integrate(self.u*self.u*self.w0))
 
    def getPrimaryElectricPotentials(self, sigma):
        """
        return the primary electric for the injections (A,B)
        """
        primary_potential={}
        self.pde.setValue(A=sigma*kronecker(3), y_dirac=Data(), Y=Data(), X=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            primary_potential[A]=self.pde.getSolution()
            lslogger.debug("primary potential for injection at %d -> %s"%(A,str(primary_potential[A])))
        return primary_potential

    def getSecondaryElectricPotentials(self, sigma, sigma0, primary_potentials):
        secondary_potentials={}
        print("getSecondaryElectricPotentials:",str(sigma))
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
        return secondary_potentials
    
    def getSigma(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m*self.alpha0
        else:
            p=m
        return self.sigma0*exp(p)
    
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)

    def _getArguments(self, m):
        
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)

        sigma=self.getSigma(ppi, isSmoothed=True)
        secondary_potentials=self.getSecondaryElectricPotentials(sigma, self.sigma0, self.phi_p)

        return secondary_potentials, sigma, ppi           
    
    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        ppi=args[2]

        mi=interpolate(m, Function(self.pde.getDomain()))
        mean=integrate(mi*self.u*self.w0)
        
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate((mi-mean*self.u)**2)
        misfit=None
        for A,B in self.data.injectionIterator():
                if misfit is None:
                    misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], self.misfit.w[(A,B)].getFunctionSpace())
                if self.useLogDefect:
                    EI_AB=length(E_AB)+self.datatol
                    misfit+=self.misfit.w[(A,B)]*log(EI_AB/self.misfit.data[(A,B)])**2
                else:
                    EI_AB=length(E_AB)
                    misfit+=self.misfit.w[(A,B)] * (1-(EI_AB/self.misfit.data[(A,B)]))**2
        A2=integrate(misfit)
        lslogger.info("sigma = %s"%(str(sigma)))
        lslogger.debug("p = %s"%(str(ppi)))
        lslogger.debug("m = %s"%(str(m)))
        lslogger.debug("L2, H1, misfit, mean= %e, %e, %e, %e"%(A0/2, A1/2, A2/2, integrate(mi)/self.VOL ))
        return (A0+A1+A2)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        
        mi=interpolate(m, Function(self.pde.getDomain()))
        mean=integrate(mi*self.u*self.w0)
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*(mi-mean*self.u)
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            EI_AB=length(E_AB)+self.datatol
            if self.useLogDefect:       
                self.pde.setValue(X=self.misfit.w[(A,B)]/(EI_AB**2)*log(EI_AB/self.misfit.data[(A,B)])*E_AB )
            else:
                raise NotImplementedError
            ustar=self.pde.getSolution()
            Y2-=inner(grad(ustar, E_AB.getFunctionSpace()),E_AB)

        if self.Spde:
            self.Spde.setValue(Y=Y2*sigma)
            Y+=self.Spde.getSolution()
        else:
            Y+=Y2*sigma/self.alpha0
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        print("XXX")
        self.Hpde.setValue(X=r[1], Y=r[0])
        #saveVTK("test", PPP=self.Hpde.getRightHandSide())
        p=self.Hpde.getSolution()
        #p*=1./Lsup(p)
        a=(1.-self.penalty)/self.penalty*integrate(self.w0*self.u* p)
        p+=a*self.u
        lslogger.debug("inverse Hessian called. search direction = %s (a= %s)"%(p,a))
        return p

class FieldInversion(MeteredCostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., sigma0=.001, region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, useLogDefect=True):
        """
        cost funtion for ERT inversion 
        
        :domain: pde domain
        :data: data, is ERTSurveyData object supporting makePrediction
        :w0: weighting L2 regularization
        :w1: weighting H1 regularization
        :sigma0: reference conductivity
        :region_fixed: mask for fixed conductivities
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        """
        super(FieldInversion, self).__init__()
        self.datatol=1e-30
        self.sigma0=sigma0
        self.stationsFMT=stationsFMT
        self.useLogDefect=useLogDefect

        if self.useLogDefect:
            lslogger.info("Misfit is using logarithm.") 
        else:
            lslogger.info("Misfit is using norm relative difference.")
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.data=data
        
        
        # when points are adjusted:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        lslogger.info("building misfit weighting (this will take some time)")
        Xcenter=ReducedFunction(domain).getX()
        X=Function(domain).getX()
        
        self.misfit = lambda: None
        self.misfit.data={}
        self.misfit.w={}

        for AB in self.data.injectionIterator(): 
            self.misfit.w[AB]=Scalar(0.,X.getFunctionSpace())
            self.misfit.data[AB]=Vector(0.,X.getFunctionSpace()) 
            
        for M in self.data.getObservationElectrodes():
            L_ABS=L_stations
            if adjustStationLocationsToElementCenter:
                xs=XStations[data.getStationNumber(M)]
                adjustmax=max(adjustmax, length(xs-self.data.getStationLocation(M)))
            else:
                xs=self.data.getStationLocation(M)
            mask=whereNegative(interpolate(length(Xcenter-xs)-L_ABS, X.getFunctionSpace()))   
            for A,B in self.data.getInjections(M):
                E0, E1, E2=self.data.getFieldData((A,B,M))
                n=E0**2+E1**2+E2**2
                if n > 0:
                    self.misfit.w[(A,B)].copyWithMask(Scalar(1./n,self.misfit.w[AB].getFunctionSpace()), mask)   # 1 where data are measured @M
                    self.misfit.data[(A,B)].copyWithMask(Vector((E0, E1, E2), self.misfit.w[AB].getFunctionSpace()), mask*[1,1,1]) # data inserted  @ M 
                    #self.misfit.data[(A,B)]=self.misfit.data[(A,B)]*(1-mask)+mask*Vector((E0, E1, E2), self.misfit.w[AB].getFunctionSpace())

        lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            s=integrate(self.misfit.w[AB]*length(self.misfit.data[AB])**2)
            #print(AB, s, integrate(length(self.misfit.w[(A,B)]*self.misfit.data[AB])**2))
            #self.misfit.data[AB]+=(1-wherePositive(self.misfit.w[AB])) # one inserted to avoid division by zero in misfit            
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))

        # primary potentials:
        lslogger.info("maximal station adjustment is %e"%adjustmax)
        lslogger.info("building primary electric fields (this will take some time)")
        self.phi_p=self.getPrimaryElectricPotentials(sigma0)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain, poisson=(abs(w1)>0) ) 
        self.Hpde.setValue(A=w1*kronecker(3), D=w0, q=region_fixed)
        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain)        
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
            if not self.alpha0 > 0: 
                self.Spde.setValue(q=region_fixed)
        else:
            self.Spde=None
 
    def getPrimaryElectricPotentials(self, sigma):
        """
        return the primary electric for the injections (A,B)
        """
        primary_potential={}
        self.pde.setValue(A=sigma*kronecker(3), y_dirac=Data(), Y=Data(), X=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            primary_potential[A]=self.pde.getSolution()
            lslogger.debug("primary potential for injection at %d -> %s"%(A,str(primary_potential[A])))
        return primary_potential

    def getSecondaryElectricPotentials(self, sigma, sigma0, primary_potentials):
        secondary_potentials={}
        print("getSecondaryElectricPotentials: sigma=",str(sigma))
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potentials[A]=self.pde.getSolution()
        return secondary_potentials
    
    def getSigma(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m*self.alpha0
        else:
            p=m
        return self.sigma0*exp(p)
    
    def _getDualProduct(self, m, r):
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        return Lsup(m)

    def _getArguments(self, m):
        
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-LOGCLIP, maxval=LOGCLIP)

        sigma=self.getSigma(ppi, isSmoothed=True)
        secondary_potentials=self.getSecondaryElectricPotentials(sigma, self.sigma0, self.phi_p)

        return secondary_potentials, sigma, ppi           
    
    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        ppi=args[2]

        mi=interpolate(m, Function(self.pde.getDomain()))
        
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit=None
        for A,B in self.data.injectionIterator():
                if misfit is None:
                    misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], self.misfit.w[(A,B)].getFunctionSpace())
                diff=self.misfit.data[(A,B)]-E_AB
                misfit+=self.misfit.w[(A,B)] * length(diff)**2
        A2=integrate(misfit)
        lslogger.info("sigma = %s"%(str(sigma)))
        lslogger.debug("p = %s"%(str(ppi)))
        lslogger.debug("m = %s"%(str(m)))
        lslogger.debug("L2, H1, misfit= %e, %e, %e"%(A0/2, A1/2, A2/2))
        return (A0+A1+A2)/2

    def _getGradient(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potentials=args[0]
        sigma=args[1]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*interpolate(m, X.getFunctionSpace())
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB=-grad(secondary_potentials[A]-secondary_potentials[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            diff=self.misfit.data[(A,B)]-E_AB
            self.pde.setValue(X=self.misfit.w[(A,B)]*diff )
            ustar=self.pde.getSolution()
            Y2+=inner(grad(ustar, E_AB.getFunctionSpace()),E_AB)

        if self.Spde:
            self.Spde.setValue(Y=Y2*sigma)
            Y+=self.Spde.getSolution()
        else:
            Y+=Y2*sigma/self.alpha0
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        self.Hpde.setValue(X=r[1], Y=r[0])
        #saveVTK("test", PPP=self.Hpde.getRightHandSide())
        p=self.Hpde.getSolution()
        #p*=1./Lsup(p)
        lslogger.debug("inverse Hessian called. search direction = %s",p)
        return p
