"""
tools for ERT inversion including ERT cost functions for inversion

by l.gross@uq.edu.au, 2021
"""

from esys.escript import *
import numpy as np
from esys.escript.minimizer import CostFunction, MinimizerException
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag, getInfLocator
import logging
from esys.weipa import saveVTK, saveSilo
from .tools import setupERTPDE


lslogger=logging.getLogger('inv.minimizer')

class DCInversionByFieldIntensity(CostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0., sigma0=.001, region_fixed=Data(), stationsFMT="e%s", 
                 weightLogDefect=0.5, adjustStationLocationsToElementCenter=True, logclip=15):
        """
        cost function for electric field intensity inversion. 
        regularization is 
        
        int( w0* m^2 + w1*grad(m)^2) where log(sigma/sigma0)=p is given a alpha0*p+alpha1*laplace(p)=m
        
        :domain: pde domain
        :data: data, is `fingal.SurveyData`, requires 'E' and - if available - 'RELERR_E' 
        :sigma0: reference conductivity
        :w0: weighting L2 regularization  int m^2
        :w1: weighting H1 regularization  int grad(m)^2
        :alpha0: regularization factor 
        :alpha1: regularization factor 
        :weightLogDefect: weighting factor for the logarithm defect in the cost funtion.
        :region_fixed: mask for fixed conductivity. needs to be set if w1>0 and w0=0 or alpha1> and alpha0=0
        :adjustStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        :L_stations: radius of electric field averaging. 
        """
        super(DCInversionByFieldIntensity, self).__init__()
        assert weightLogDefect >=0 and weightLogDefect<=1, "weightLogDefect needs to be between 0 and 1."
        self.datatol=1e-30
        self.sigma0=sigma0
        self.stationsFMT=stationsFMT
        self.weightLogDefect=weightLogDefect
        self.logclip=logclip

        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        self.data=data


        # when points are adjusted get the element center locations:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        if getMPIRankWorld() == 0: lslogger.info("building misfit weighting (this will take some time)")
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
                E=self.data.getFieldIntensityData((A,B,M))
                RELERR=self.data.getFieldIntensityRelError((A,B,M))
                self.misfit.w[(A,B)].copyWithMask(Scalar(1/RELERR**2,self.misfit.w[AB].getFunctionSpace()), mask)   # >0 where data are measured @M
                self.misfit.data[(A,B)].copyWithMask(Scalar(E,self.misfit.w[AB].getFunctionSpace()), mask) # data inserted  @ M 

        if getMPIRankWorld() == 0: lslogger.debug("re-scaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            s=integrate(self.misfit.w[AB])
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # data inserted  @ M 
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        if getMPIRankWorld() == 0: 
            lslogger.info("maximal station adjustment is %e"%adjustmax)
            lslogger.info("building primary electric fields (this will take some time)")
        
        self.phi_p=self.getPrimaryElectricPotentials(sigma0)
        lslogger.info("Primary potentials for %d injections calculated."%(len(self.phi_p) ))
        # this defines the regularization:
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain) 
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
        return the primary electric potential for all injections (A,B) using conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_potential
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
            txt=str(primary_potential[A])
            if getMPIRankWorld() == 0: 
                lslogger.debug("primary potential for injection at %d -> %s"%(A,txt))
        return primary_potential

    def getSecondaryElectricPotentials(self, sigma, sigma0, primary_potentials):
        """
        return the primary electric  potential for all injections (A,B) using conductivity sigma
        for the primary conductivity sigma0 and potentials 
        :sigma: (primary) conductivity distribution primary_potentials
        :return: dictonary of injections (A,B)->secondary_potential
        """
        
        secondary_potential={}
        txt=str(sigma)
        if getMPIRankWorld() == 0: 
            lslogger.debug("getSecondaryElectricPotentials: sigma="+txt)
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potential[A]=self.pde.getSolution()
        return secondary_potential
    
    def optimizeSigma0(self, m):
        """
        returns a new conductivity, a scaling factor and new misfit by minimizing the misfit using conductivity sigma=f*sigma0 over 
        factor f.
        """
        raise NotImplemented
    def scaleSigma0(self, f=1.):
        """
        rescales sigma0 by factor f. 
        """
        raise NotImplemented
        
    def getSigma(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m. If isSmoothed=True 
        it is assumed that m is already smoothed by (alpha0*I+alpha1*laplace)^{-1}. Otherwise smoothing is applied.
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
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """        
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """        
        return Lsup(m)

    def _getArguments(self, m):
        """
        returns values that are used for both the forward as well as the gradient calculation
        """        
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)

        sigma=self.getSigma(ppi, isSmoothed=True)
        
        txt1, txt2=str(ppi), str(m)
        if getMPIRankWorld() == 0: 
            lslogger.debug("p = %s"%(txt1))
            lslogger.debug("m = %s"%(txt2))

        secondary_potential=self.getSecondaryElectricPotentials(sigma, self.sigma0, self.phi_p)
        return secondary_potential, sigma, ppi           
    
    def _getValue(self, m, *args):
        """
        return the value of the cost function. Overwrites `getValue` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        
        secondary_potential=args[0]
        sigma=args[1]
        ppi=args[2]

        mi=interpolate(m, Function(self.pde.getDomain()))
        
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit_log=None
        for A,B in self.data.injectionIterator():
                if misfit_log is None:
                    misfit_log=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                    misfit_quad=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

                E_AB=-grad(secondary_potential[A]-secondary_potential[B]+self.phi_p[A]-self.phi_p[B], self.misfit.w[(A,B)].getFunctionSpace())
                EI_AB=length(E_AB)+self.datatol
                
                misfit_log+=self.misfit.w[(A,B)]*(log(EI_AB/self.misfit.data[(A,B)]))**2
                misfit_quad+=self.misfit.w[(A,B)] * (1-(EI_AB/self.misfit.data[(A,B)]))**2

        A2=integrate(misfit_log)
        A3=integrate(misfit_quad)
        if getMPIRankWorld() == 0: 
            lslogger.debug("L2, H1, misfit quad, log= %e, %e, %e, %e"%(A0/2, A1/2, A3/2, A2/2))
        return (A0+A1+(1-self.weightLogDefect)*A3+self.weightLogDefect*A2)/2

    def _getGradient(self, m, *args):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potential=args[0]
        sigma=args[1]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*interpolate(m, X.getFunctionSpace())
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB=-grad(secondary_potential[A]-secondary_potential[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            EI_AB=length(E_AB)+self.datatol
            D=self.misfit.data[(A,B)]
            m_log=log(EI_AB/D)
            m_quad=1-(EI_AB/D)
    
            self.pde.setValue(X=(self.misfit.w[(A,B)]*(m_log/(EI_AB**2)*self.weightLogDefect - m_quad/(D*EI_AB)*(1-self.weightLogDefect)) ) * E_AB)
            ustar=self.pde.getSolution()
            Y2-=inner(grad(ustar, E_AB.getFunctionSpace()),E_AB)

        if self.Spde:
            self.Spde.setValue(Y=Y2*sigma)
            Y+=self.Spde.getSolution()
        else:
            Y+=Y2*sigma/self.alpha0
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        txt=str(p)
        if getMPIRankWorld() == 0: 
            lslogger.debug("inverse Hessian called. search direction = %s",txt)
        return p
    
class ChargeabilityInversionByField(CostFunction):
    """
    cost function for electric field intensity  inversion (aka FullWaver) 
    """
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, L_stations=1., w0=0., w1=1., alpha0=1., alpha1=0.,  gamma0=0.001, sigma=0.001, region_fixed=Data(), stationsFMT="e%s", 
                 adjustStationLocationsToElementCenter=True, logclip=15, weightLogDefect=0.):
        """
        cost function for chargeability inversion based electric fields. 
        regularization is 
        
        int( w0* m^2 + w1*grad(m)^2) where log(sigma/sigma0)=p is given a alpha0*p+alpha1*laplace(p)=m
        
        :domain: pde domain
        :data: data, is `fingal.SurveyData`, requires 'GAMMA' if available  'RELERR_GAMMA' 
        :gamma0: reference modified chargeability
        :sigma: conductivity
        :w0: weighting L2 regularization  int m^2
        :w1: weighting H1 regularization  int grad(m)^2
        :alpha0: regularization factor 
        :alpha1: regularization factor 
        :weightLogDefect: weighting factor for the logarithm defect in the cost funtion.
        :region_fixed: mask for fixed conductivity. needs to be set if w1>0 and w0=0 or alpha1> and alpha0=0
        :adjustStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format used to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        :L_stations: radius of electric field averaging. 
        """
        super(ChargeabilityInversionByField, self).__init__()
        self.datatol=1e-30
        self.logclip=logclip
        if getMPIRankWorld() == 0:
             if weightLogDefect>0:
                lslogger.info("weightLogDefect>0 but ignored.")                 
             lslogger.info("building misfit weighting (this will take some time)")

             
        self.weightLogDefect=weightLogDefect
        self.sigma=sigma
        self.gamma0=gamma0
        self.useDifferenceOfFields=False
        self.stationsFMT=stationsFMT
        self.misfitFunctionSpace=Function(domain)
        # setup PDE:
        self.pde=setupERTPDE(domain)
        
        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

         
        self.data=data
        
        
        # when points are adjusted to match element centers:
        adjustmax=0.
        if adjustStationLocationsToElementCenter:
            station_locations=[]
            for s in data.getStationNumeration():
                station_locations.append(data.getStationLocation(s))
            XStations=Locator(ReducedFunction(domain), station_locations).getX()

        ########################################################
        Xcenter=ReducedFunction(domain).getX()
        X=self.misfitFunctionSpace.getX()
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
                GAMMA=self.data.getModifiedChargeabilityData((A,B,M))
                RELERR=self.data.getModifiedChargeabilityRelError((A,B,M))
                if abs(GAMMA) > 0.:
                    self.misfit.w[(A,B)].copyWithMask(Scalar(1./RELERR**2,self.misfitFunctionSpace), mask)   # 1 where data are measured @M
                    self.misfit.data[(A,B)].copyWithMask(Scalar(GAMMA,self.misfitFunctionSpace), mask) # data inserted  @ M 
                    
        if adjustStationLocationsToElementCenter and getMPIRankWorld() == 0:
            lslogger.info("maximal station adjustment is %e"%adjustmax)
        if getMPIRankWorld() == 0:
            lslogger.debug("rescaling of misfit weights:")
        for AB in self.data.injectionIterator(): 
            self.misfit.data[AB]+=whereNonPositive(self.misfit.w[AB]) # insert 1's to avoid division by zero
            s=integrate(self.misfit.w[AB])
            assert s>0, "no observation for dipole %s. Maybe you need to increase the value for L_stations."%(str(AB))
            if s > 0:
                self.misfit.w[AB]*=1./(s*len(self.misfit.w))
        # primary potentials:
        self.phi_p=self.getElectricPotentials(self.sigma)
        #self.secondary_potential=self.getSecondaryElectricPotentials(self.sigma, self.sigma0, self.phi_p)
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
 
    def getElectricPotentials(self, sigma):
        """
        return the primary electric potentials for the injections (A,B) for conductivity sigma
        """
        potential={}
        self.pde.setValue(A=sigma*kronecker(3), X=Data(), Y=Data(), y_dirac=Data())
        for A in self.data.getListOfInjectionStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(A,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%A,1.)
            self.pde.setValue(y_dirac=s)
            potential[A]=self.pde.getSolution()
            txt=str(potential[A])
            if getMPIRankWorld() == 0: 
                lslogger.debug("primary potential for injection at %d -> %s"%(A,txt))
        if getMPIRankWorld() == 0:
            lslogger.debug("primary potential for %s injection calculated."%len(potential))
        return potential

    def getSecondaryElectricPotentials(self, gamma, primary_potentials):
        """
        return the secondary electric potentials for the injections (A,B) for conductivity sigma/(1+gamma)
        """        
        secondary_potential={}
        self.pde.setValue(A=self.sigma/(1+gamma)*kronecker(3), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=self.sigma*gamma/(1+gamma)*grad(primary_potentials[A]))
            secondary_potential[A]=self.pde.getSolution()
        if getMPIRankWorld() == 0:
            lslogger.debug("%s secondary potential calculated"%(len(secondary_potential)))
        return secondary_potential
            
    def getChargeability(self, m, isSmoothed=False):
        """
        return chargeability (eta) for a given property function m.
        If isSmoothed=True it is assumed that m is already smoothed by (alpha0*I+alpha1*laplace)^{-1}. Otherwise smoothing is applied.
        """
        gamma=self.getGamma(m, isSmoothed)
        return gamma/(1.+gamma)
    
    def getGamma(self, m, isSmoothed=False):
        """
        return modified chargeability (gamma) for a given property function m
        If isSmoothed=True it is assumed that m is already smoothed by (alpha0*I+alpha1*laplace)^{-1}. Otherwise smoothing is applied.
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m/self.alpha0
        else:
            p=m
        gamma=self.gamma0*exp(clip(p, minval=-self.logclip, maxval=self.logclip))
        return gamma
        
    def _getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """    
        return integrate(r[0]*m + inner(r[1], grad(m)))
    
    def _getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """    
        return Lsup(m)


    def _getArguments(self, m):
        """
        returns values that are used for both the forward as well as the gradient calculation
        """   
        if self.Spde :
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
        else:
            p=m/self.alpha0
                
        gamma=self.getGamma(p, isSmoothed=True)
        
        gammai=interpolate(gamma,  Function(self.pde.getDomain()))
        secondary_potential=self.getSecondaryElectricPotentials(gammai, self.phi_p)
        return gammai, secondary_potential, 

    def _getValue(self, m, *args):
        """
        return the value of the cost function. Overwrites `getValue` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        gammai=args[0]
        secondary_potential=args[1]

        mi=interpolate(m,  Function(self.pde.getDomain()))
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        
        misfit=Scalar(0.,self.misfitFunctionSpace )
        for A, B in self.data.injectionIterator():
            E =-grad(self.phi_p[A]-self.phi_p[B], misfit.getFunctionSpace())
            DE=-grad(secondary_potential[A]-secondary_potential[B], misfit.getFunctionSpace())
            L_E2=length(E)**2
            mfquad=1 - safeDiv( inner(DE, E), self.misfit.data[(A,B)]*L_E2)
            misfit+=self.misfit.w[(A,B)]*mfquad**2
        A2=integrate(misfit)

        if lslogger.isEnabledFor(logging.DEBUG):
            strgamma=str(gammai)
            strm=str(mi)
            if getMPIRankWorld() == 0:
                lslogger.debug("gamma = %s"%strgamma)
                lslogger.debug("m = %s"%strm)
                lslogger.debug("L2, H1, misfit quad = %e, %e, %e"%(A0/2, A1/2, A2/2))

        return (A0+A1+A2)/2

    def _getGradient(self, m, *args):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        gammai=args[0]
        secondary_potential=args[1]
        mi=interpolate(m,  Function(self.pde.getDomain()))

        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*mi
        self.pde.setValue(A=self.sigma/(1+gammai)*kronecker(3), X=Data(), Y=Data(), y_dirac=Data())
        Y2=Scalar(0.,self.misfitFunctionSpace )
        for A, B in self.data.injectionIterator():
            E =-grad(self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
            DE=-grad(secondary_potential[A]-secondary_potential[B], Y2.getFunctionSpace())
            L_E2=length(E)**2
            mfquad=1 - safeDiv( inner(DE, E), self.misfit.data[(A,B)]*L_E2)

            #self.pde.setValue(X=self.misfit.w[(A,B)]*(self.misfit.data[(A,B)]*L_E2 - inner(DE, E))/(L_E2+self.datatol**2)**2* E)
            self.pde.setValue(X=self.misfit.w[(A,B)]*safeDiv(mfquad, self.misfit.data[(A,B)]*L_E2) * E)
            ustar=self.pde.getSolution()
            Y2+=-inner(grad(ustar, E.getFunctionSpace()),E+DE)

        if self.Spde:
            self.Spde.setValue(Y=Y2*gammai/(1+gammai)**2*self.sigma)
            Y+=self.Spde.getSolution()
        else:                
            Y+=Y2*gammai*self.sigma/self.alpha0

        
        return ArithmeticTuple(Y, X)

        
    def _getInverseHessianApproximation(self, m, r, *args):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        txt=str(p)
        if getMPIRankWorld() == 0: 
            lslogger.debug("inverse Hessian called. search direction = %s",txt)
        return p


class DCInversionByField(CostFunction):
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
        secondary_potential={}
        print("getSecondaryElectricPotentials: sigma=",str(sigma))
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for A in primary_potentials:
            self.pde.setValue(X=(sigma0-sigma)*grad(primary_potentials[A]))
            secondary_potential[A]=self.pde.getSolution()
        return secondary_potential
    
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
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)

        sigma=self.getSigma(ppi, isSmoothed=True)
        secondary_potential=self.getSecondaryElectricPotentials(sigma, self.sigma0, self.phi_p)

        return secondary_potential, sigma, ppi           
    
    def _getValue(self, m, *args):
        if len(args)==0:
            args=self.getArguments(m)
        secondary_potential=args[0]
        sigma=args[1]
        ppi=args[2]

        mi=interpolate(m, Function(self.pde.getDomain()))
        
        A1=self.w1*integrate(length(grad(m))**2)
        A0=self.w0*integrate(mi**2)
        misfit=None
        for A,B in self.data.injectionIterator():
                if misfit is None:
                    misfit=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )
                E_AB=-grad(secondary_potential[A]-secondary_potential[B]+self.phi_p[A]-self.phi_p[B], self.misfit.w[(A,B)].getFunctionSpace())
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
        secondary_potential=args[0]
        sigma=args[1]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*interpolate(m, X.getFunctionSpace())
        self.pde.setValue(A=sigma*kronecker(self.pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        Y2=None
        for A, B in self.data.injectionIterator():
            if Y2 is None:
                Y2=Scalar(0.,self.misfit.w[(A,B)].getFunctionSpace() )

            E_AB=-grad(secondary_potential[A]-secondary_potential[B]+self.phi_p[A]-self.phi_p[B], Y2.getFunctionSpace())
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
