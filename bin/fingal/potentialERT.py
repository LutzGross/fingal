"""
tools for ERT inversion including ERT cost functions for inversion

by l.gross@uq.edu.au, 2021
"""

from esys.escript import *
from esys.downunder import MeteredCostFunction
from .tools import setupERTPDE

import numpy as np
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import Locator, ArithmeticTuple, MaskFromTag, getInfLocator
import logging
from esys.weipa import saveVTK, saveSilo

lslogger=logging.getLogger('inv.minimizer')

class PotentialERT(MeteredCostFunction):
    provides_inverse_Hessian_approximation=True
    def __init__(self, domain, data, w0=0., w1=1., alpha0=1., alpha1=0., sigma0=.01, region_fixed=Data(), stationsFMT="e%s", weightLogDefect=0., logclip=15):
        """
        cost function for ERT inversion based on potential data. The unknown of the inversion is
        
        m=alpha0*p - alpha1 div(grad(p)) with p=log(sigma/sigma0)
      
        where sigma is the unknown conductivity and sigma0 an assumed conductivity.
        
        the regularization R(m) is a combination of L2 and H1:
        
        R(m) = 1/2 * integrate( w0*m**2 + w1* |grad(m)|^2 )
        
        The misfit is given quadratic 
        
        :param domain: inversion domain
        :param data: data, is ERTSurveyData object supporting makeResistencePrediction and getResistenceData
        :param w0: weighting for L2 regularization
        :param w1: weighting for H1 regularization
        :param sigma0: reference conductivity
        :param region_fixed: mask for fixed conductivity distribution (this where the value sigma0 is used)
        :param stationsFMT: format used to map station keys k to mesh tags stationsFMT%k 
        :param logclip: values of log(sigma/sigma0)>logclip are set to logclip to avoid overflow
        """
        super(PotentialERT, self).__init__()

        assert weightLogDefect <= 1 and weightLogDefect >=0, "weightLogDefect  needs to be between 0 and 1"
        self.stationsFMT=stationsFMT
        self.w0=w0
        self.w1=w1
        self.alpha0=alpha0
        self.alpha1=alpha1
        self.setSigma0(sigma0)
        self.data=data
        self.logclip=logclip
        self.sfactor=1./float(self.data.getNumObservations())
        self.weightLogDefect=weightLogDefect
        
        self.pde=setupERTPDE(domain)

        x=self.pde.getDomain().getX()[0]
        y=self.pde.getDomain().getX()[1]
        z=self.pde.getDomain().getX()[2]
        self.pde.setValue(q=whereZero(x-inf(x))+whereZero(x-sup(x))+ whereZero(y-inf(y))+whereZero(y-sup(y))+whereZero(z-inf(z)))

        # used for Hessian inverse
        self.Hpde=setupERTPDE(domain)
        self.Hpde.setValue(A=self.w1*kronecker(3), D=self.w0, q=region_fixed)

        if self.alpha1 > 0: 
            self.Spde=setupERTPDE(domain)        
            self.Spde.setValue(A=self.alpha1*kronecker(3), D=self.alpha0)
        else:
            self.Spde=None
            
 
        station_locations=[]
        for s in self.data.getStationNumeration():
            station_locations.append(self.data.getStationLocation(s))
        self.locators=Locator(Solution(domain), station_locations)
   
        
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
        # smooth m by solving m = alpha0*p - alpha1 div(grad(p)) 
        if self.Spde:
            self.Spde.setValue(Y=m)
            p=self.Spde.getSolution()
            ppi=clip(interpolate(p, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)
        else:
            ppi=clip(interpolate(m/self.alpha0, Function(self.pde.getDomain())), minval=-self.logclip, maxval=self.logclip)
        # now calculate the new sigma
        sigma=self.getSigma(ppi, isSmoothed=True)
        if getMPIRankWorld() == 0:
            lslogger.info("sigma = %s"%str(sigma))
        
        # now we solve for all electrodes creating a dictionary `potentials[ip]` and potentials_at_stations[ip] where `ip` is the electrode id.
        # note that we not distinguish between charging and measuring electrodes.
        self.pde.setValue(A=sigma*kronecker(3))
        potentials={}
        potentials_at_stations={}
        for ip in self.data.getListOfAllStations():
            s=Scalar(0.,DiracDeltaFunctions(self.pde.getDomain()))
            if self.stationsFMT is None:
                s.setTaggedValue(ip,1.)
            else:    
                s.setTaggedValue(self.stationsFMT%ip,1.)
            self.pde.setValue(y_dirac=s)
            u=self.pde.getSolution()
            potentials[ip]=u
            potentials_at_stations[ip]=np.array(self.locators(u))
        if getMPIRankWorld() == 0:
            lslogger.info("Potential for %d electrodes calculated."%len(potentials))
        dV=self.data.makeResistencePrediction(values=potentials_at_stations)

        return sigma, potentials, dV
    
    def optimizeSigma0(self, m=0):
        """
        this returns a corrected value (float) for sigma0 that gives a better initial data match 
        """
        sigma, potentials, dV=self.getArguments(m)
        A1=0.
        A2=0.
        for t in self.data.tokenIterator(): 
            v=dV[t]
            d=self.data.getResistenceData(t)
            e=self.data.getResistenceRelError(t)
            if abs(d) > 0:
                A1+=v/(d*e**2)
                A2+=(v/(d*e))**2
        if A2 > 0 and A1 >0:
            f_opt=A1/A2    
        else:
            f_opt=1.
        sigma_opt=self.sigma0/f_opt

        defect=0.
        for t in self.data.tokenIterator(): # t=(A,B,M,N) (or so)
            v=dV[t]
            d=self.data.getResistenceData(t)
            e=self.data.getResistenceRelError(t)
            defect+=((f_opt*v-d)/(e*d))**2
        defect*=self.sfactor
        return sigma_opt, 1./f_opt, defect
    
    def setSigma0(self, sigma0=0.01):
        """
        set a new value for sigma0
        """
        self.sigma0=sigma0
        
    def getSigma(self, m, isSmoothed=False):
        """
        return the conductivity for a given property function m
        """
        if not isSmoothed:
            if self.Spde :
                self.Spde.setValue(Y=m)
                p=self.Spde.getSolution()
            else:
                p=m/self.alpha0
        else:
            p=m
        return self.sigma0*exp(p)
    
    def _getValue(self, m, *args):
        """
        return the value of the cost function. Overwrites `getValue` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        sigma=args[0]    
        potentials=args[1]
        dV=args[2]
        
        # regularization terms:
        A1=integrate(self.w1*length(grad(m))**2)
        A0=integrate(self.w0*interpolate(m, Function(self.pde.getDomain()))**2)

        # misfit
        A2, A3=0., 0.
        for t in self.data.tokenIterator(): # t=(A,B,M,N) (or so)
            v=dV[t]
            d=self.data.getResistenceData(t)
            e=self.data.getResistenceRelError(t)
            difflog=log(abs(v/d))/e
            diffquad=(v-d)/(e*d)
                
            A2+=difflog**2
            A3+=diffquad**2

        A2*=self.sfactor
        A3*=self.sfactor
        if getMPIRankWorld() == 0:
            lslogger.info("reg: L0, H1, misfit: log, quad = %e, %e, %e, %e"%(A0/2, A1/2, A2/2, A3/2))
        
        return (A0+A1+(1-self.weightLogDefect)*A3+self.weightLogDefect*A2)/2
 

    def _getGradient(self, m, *args):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        if len(args)==0:
            args=self.getArguments(m)
        sigma=args[0]    
        potentials=args[1]
        dV=args[2]
        
        # gradient of the regularization part:
        X=self.w1*grad(m)
        Y=self.w0*interpolate(m,X.getFunctionSpace())

        #defects={}
        #for s in self.data.injectionIterator(): # s=(A,B)
            #defects[s]=Scalar(0, DiracDeltaFunctions(self.pde.getDomain()))
        Y2=Scalar(0., Y.getFunctionSpace())
        for inj in self.data.injectionIterator():
            if self.data.hasDipoleInjections():
                u=potentials[inj[0]]-potentials[inj[1]]
                idx=2
            else:
                u=potentials[inj]
                idx=1
            ustar=Scalar(0, u.getFunctionSpace())
            for t in self.data.getObservations(inj, insertSource=True):
                    v=dV[t]
                    d=self.data.getResistenceData(t)
                    e=self.data.getResistenceRelError(t)                    

                    difflog=safeDiv( log(abs(v/d)) , v*e**2)
                    diffquad=(v-d)/(e*d)**2                     
                    diff=((1-self.weightLogDefect)*diffquad+self.weightLogDefect*difflog)*self.sfactor
                    if self.data.hasDipoleMeasurements():
                        M,N=t[idx:]
                        ustar+=(potentials[M]-potentials[N])*diff
                    else:
                        M=t[idx]
                        ustar+=potentials[M]*diff
            Y2+=inner(grad(ustar),grad(u))
        Y2*=-sigma

        if self.Spde:
            self.Spde.setValue(Y=Y2)
            p=self.Spde.getSolution()
            Y+=interpolate(p, Function(self.pde.getDomain()))
        else:
            Y+=Y2/self.alpha0
            
        return ArithmeticTuple(Y, X)
        
    def _getInverseHessianApproximation(self, m, r, *args):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        p=self.Hpde.getSolution()
        return p
