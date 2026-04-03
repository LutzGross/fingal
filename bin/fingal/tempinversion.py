#
#
#
from esys.escript import *
from esys.weipa import saveSilo
from esys.escript.minimizer import CostFunction, MinimizerException

from .tools import setupERTPDE
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple
from .inversionsIP import IPMisfitCostFunction

class IPConductivityModelTemplate(object):
    """
    this a template for class for providing an electric conductivity model for IP inversion
    the basic functionality is that it provides values for secondary conductivity sigma_0 as function of normalized chargeability M_n but may extend
    to other models such as temperature dependencies
    """

    def __init__(self, T_min=1, T_max=105000, **kargs):
        """
        :param T_min: lower cut-off temperature
        :param T_max: upper cut-off temperature

        """
        self.T_min = T_min
        self.T_max = T_max

    def getDCConductivity(self, Mn):
        raise NotImplementedError

    def getDsigma_0DMn(self, Mn):
        raise NotImplementedError

    def getDMnDsigma_0(self, sigm0):
        raise NotImplementedError


class ConductivityModelByTemperature(IPConductivityModelTemplate):
    """
    temperature based conductivity model
    """

    def __init__(self, sigma_0_ref=1.52e-4, Mn_ref=5.054e-5, T_ref=15, P_0=2.9, P_Mn=2.5, T_min=1, T_max=105000, **kargs):
        """
        :param T_ref: reference temperature
        :param sigma_0_ref: secondary conductivity at reference temperature
        :param Mn_ref: chargeability at at reference temperature
        :param P_0: exponent of secondary conductivity  temperature model
        :param P_Mn: exponent of chargeability temperature model
        :param T_min: lower cut-off temperature
        :param T_max: upper cut-off temperature
        """
        super().__init__(T_min=T_min, T_max=T_max, **kargs)
        self.T_ref = T_ref
        self.sigma_0_ref = sigma_0_ref
        self.Mn_ref = Mn_ref
        self.P_0 = float(P_0)
        self.P_Mn = float(P_Mn)

    def getInstantanousConductivity(self, T):
        """
        return the secondary conductivity sigma_0 for given temperature T or given chargeability Mn
        """
        sigma_oo = self.getChargeability(T=T) + self.getDCConductivity(T=T)
        return sigma_oo

    def getDCConductivity(self, T=None, Mn=None):
        """
        return the secondary conductivity sigma_0 for given temperature T or given chargeability Mn
        """
        if not T is None:
            sigma_0 = self.sigma_0_ref * (T / self.T_ref) ** self.P_0
        elif not Mn is None:
            sigma_0 = self.sigma_0_ref * (Mn / self.Mn_ref) ** (self.P_0 / self.P_Mn)
        return sigma_0

    def getTemperature(self, Mn=None, sigma_0=None):
        """
        return the temperature for  given secondary conductivity sigma_0  or given  chargeability Mn
        """
        if not Mn is None:
            T = self.T_ref * (Mn / self.Mn_ref) ** (1. / self.P_Mn)
        elif not sigma_0 is None:
            T = self.T_ref * (sigma_0 / self.sigma_0_ref) ** (1. / self.P_0)
        return T

    def getChargeability(self, T=None, sigma_0=None):
        """
        return the chargeability for given secondary conductivity sigma_0  or for given temperature T
        """
        if not T is None:
            Mn = self.Mn_ref * (T / self.T_ref) ** self.P_Mn
        elif not sigma_0 is None:
            Mn = self.Mn_ref * (sigma_0 / self.sigma_0_ref) ** (self.P_Mn / self.P_0)
        return Mn

    def getDsigma_0DMn(self, Mn):
        Dsigma_0DMn = (self.P_0 / self.P_Mn) * (self.sigma_0_ref / self.Mn_ref) * (Mn / self.Mn_ref) ** (
                    self.P_0 / self.P_Mn - 1)
        return Dsigma_0DMn

    def getDMnDsigma_0(self, sigma_0):
        DMnDsigma_0 = (self.P_Mn / self.P_0) * (self.Mn_ref / self.sigma_0_ref) * (sigma_0 / self.sigma_0_ref) ** (
                    self.P_Mn / self.P_0 - 1)
        return DMnDsigma_0

    def getDsigma_0DT(self, T):
        Dsigma_0DT = self.P_0 * (self.sigma_0_ref / self.T_ref) * (T / self.T_ref) ** (self.P_0 - 1)
        return Dsigma_0DT

    def getDMnDT(self, T):
        DMnDT = self.P_Mn * (self.Mn_ref / self.T_ref) * (T / self.T_ref) ** (self.P_Mn - 1)
        return DMnDT



class InversionIPByFlux(IPMisfitCostFunction):
    """

    """

    def __init__(self, domain, data, maskZeroPotential = None,
                 conductivity_model=IPConductivityModelTemplate(),
                 surface_temperature = None, mask_surface_temperature = None,
                 sigma_src=None, pde_tol=1e-8, stationsFMT="e%s", length_scale=None,
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, w1=1, reg_tol=None,
                 thermal_conductivity=1, heat_production=None, surface_flux=None,
                 logger=None, **kargs):
        """
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.
        :param conductivity_model: conductivity model
        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
        :

        :param pde_tol: tolerance for solving the forward PDEs and recovering m from M.
        :param stationsFMT: format string to convert station id to mesh label
        :param m_ref: reference property function
        :param zero_mean_m: constrain m by zero mean.
        :param useLogMisfitDC: if set logarithm of DC data is used in misfit.
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolDC: relative tolerance for damping small DC data on misfit
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param length_scale: length scale, w0=(1/length scale)**2 is weighting factor for |grad m|^2 = |M|^2 in
                            the regularization term. If `None`, the term is dropped.
        :param theta: weigthing factor x-grad term
        :param fixTop: if set m[0]-m_ref[0] and m[1]-m_ref[1] are fixed at the top of the domain
                        rather than just on the left, right, front, back and bottom face.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param save_memory: if not set, the three PDEs for the three conponents of the inversion unknwon
                            with the different boundary conditions are held. Otherwise, boundary conditions
                            are updated for each component which requires refactorization which obviously
                            requires more time.
        :param logger: the logger, if not set, 'fingal.IPInversion.H2' is used.
        """
        if sigma_src == None:
            sigma_src = conductivity_model.sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversionByTemperature')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                         useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                         useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                         weightingMisfitDC=weightingMisfitDC,
                         logger=logger, **kargs)
        self.conductivity_model = conductivity_model
        if length_scale is None:
            self.a = None
        else:
            self.a = (1. / length_scale) ** 2
        self.logger.info("Length scale factor is %s." % (str(length_scale)))
        # PDE to recover temperature:
        if surface_temperature is None:
            surface_temperature = Scalar( self.conductivity_model.T_ref, Solution(domain))
        if mask_surface_temperature is None:
            z =Solution(domain).getX()[2]
            mask_surface_temperature = whereZero(z-sup(z))

        x = domain.getX()
        qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
        qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])))
        qz = whereZero(x[2] - inf(x[2]))

        self.Tpde = setupERTPDE(self.domain)
        self.Tpde.getSolverOptions().setTolerance(pde_tol)
        #optionsG = self.Tpde.getSolverOptions()
        #from esys.escript.linearPDEs import SolverOptions
        #optionsG.setSolverMethod(SolverOptions.DIRECT)
        self.Tpde.setValue(A = thermal_conductivity * kronecker(3),
                           q=mask_surface_temperature, r=surface_temperature)
        ## get the initial temperature:
        if not heat_production is None:
            self.Tpde.setValue( Y=heat_production)
        if not surface_flux is None:
                self.Tpde.setValue( y=surface_flux)
        self.T_background = self.Tpde.getSolution()
        self.logger.debug("Background Temperature = %s" % (str(self.T_background)))
        self.Tpde.setValue(r=Data(), Y=Data(), y=Data(), q=qx+qy)
        # ... regularization
        if not reg_tol:
            reg_tol = min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.TangentialConstraints = []

        if self.a is None:
            #self.TangentialConstraints = [qx, qy, qz]
            self.TangentialConstraints = [qy+qz, qz+qx, qx+qy]
            self.Hpde.setValue(A=kronecker(3))
        else:
            self.TangentialConstraints = [qy + qz, qz + qx, qx + qy]
            self.Hpde.setValue(A=kronecker(3), D=self.a)
        self.setW1(w1)
    #====
    def getTemperature(self, M):
        """
        returns the temperature from property function M
        """
        self.Tpde.setValue(X=M, Y =Data())
        T = self.Tpde.getSolution()
        self.logger.debug("Correction Temperature = %s" % (str(T)))
        T += self.T_background
        T = clip(T, minval=self.conductivity_model.T_min, maxval=self.conductivity_model.T_max)
        return T

    def setW1(self, w1):
        self.w1 = w1
        self.logger.debug(f'w1 = {self.w1:g}')

    def getSigma0(self, T, applyInterploation=False):
        """
        get the sigma_0 from temperature T
        """
        if applyInterploation:
            iT = interpolate(T, Function(self.domain))
        else:
            iT = T
        sigma_0 =  self.conductivity_model.getDCConductivity(T=iT)
        return sigma_0

    def getMn(self, T, applyInterploation=False):
        if applyInterploation:
            iT= interpolate(T, Function(self.domain))
        else:
            iT= T
        Mn  =  self.conductivity_model.getChargeability(T=iT)
        return Mn

    def getDsigma0DT(self, sigma_0, T):
        iT = interpolate(T, sigma_0.getFunctionSpace())
        return self.conductivity_model.getDsigma_0DT(T=iT)

    def getDMnDT(self, Mn, T):
        iT = interpolate(T, Mn.getFunctionSpace())
        return self.conductivity_model.getDMnDT(T=iT)

    def extractPropertyFunction(self, M):
        return M

    def getArguments(self, M):
        """
        precalculation
        """
        M2 = self.extractPropertyFunction(M)
        iM = interpolate(M2, Function(self.domain))
        T=self.getTemperature(iM)
        iT = interpolate(T, Function(self.domain))
        iT_stations = self.grabValuesAtStations(T)
        self.logger.debug("Flux M[0] = %s" % (str(iM[0])))
        self.logger.debug("     M[1] = %s" % (str(iM[1])))
        self.logger.debug("     M[2] = %s" % (str(iM[2])))
        self.logger.debug("Temperature = %s" % (str(iT)))
        #self.logger.debug("Temperature at stations = %s" % (str(iT_stations)))

        isigma_0 = self.getSigma0(iT)
        isigma_0_stations = self.getSigma0(iT_stations)
        iMn = self.getMn(iT)
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return iT, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, M, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)
        gM = grad(M, where=iT.getFunctionSpace())
        iM = interpolate(M, iT.getFunctionSpace())
        if self.a:
            R = 1. / 2. * integrate(self.w1 * (length(gM) ** 2 + self.a * length(iM) ** 2) )
        else:
            R = 1. / 2. * integrate(self.w1 * length(gM) ** 2  )

        V = R + misfit_DC + misfit_IP
        self.logger.debug(
                f'misfit ERT, IP; reg, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e} = {V:e}')
        self.logger.debug(
                f'ratios ERT, IP; reg \t=  {misfit_DC / V * 100:g}, {misfit_IP / V * 100:g};  {R / V * 100:g}')

        return V

    #=======================================
    def getGradient(self, M, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gM = grad(M, where=iT.getFunctionSpace())
        iM = interpolate(M, iT.getFunctionSpace())

        X = self.w1 * gM
        if self.a:
            Y = self.w1 * self.a * iM
        else:
            Y = Vector(0., iT.getFunctionSpace() )
        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)
        Dsigma_0DT = self.getDsigma0DT(isigma_0, iT)
        DMnDT = self.getDMnDT(iMn, iT)
        Ystar = DMisfitDMn * DMnDT + DMisfitDsigma_0 * Dsigma_0DT
        self.Tpde.setValue(Y = Ystar, X=Data())
        Tstar = self.Tpde.getSolution()
        Y += grad(Tstar, where=Y.getFunctionSpace())
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, M, iT, isigma_0, isigma_0_stations, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        P = Data(0., (3,), Solution(self.domain))
        if self.TangentialConstraints:
            for k in [0, 1, 2]:
                self.Hpde.setValue(X=r[1][k], Y=r[0][k], q=self.TangentialConstraints[k])
                P[k] = self.Hpde.getSolution() / self.w1
                self.logger.debug(f"search direction component {k} = {str(P[k])}.")

        else:
            for k in [0, 1, 2]:
                self.Hpde.setValue(X=r[1][k], Y=r[0][k])
                P[k] = self.Hpde.getSolution() / self.w1
                self.logger.debug(f"search direction component {k} = {str(P[k])}.")
        return P

    def getDualProduct(self, M, r):
        """
        dual product of gradient `r` with increment `V`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0], M) + inner(r[1], grad(M)))

    def getNorm(self, M):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(M)

    def getSqueezeFactor(self, M, p):
        return None

class InversionIPBySource(IPMisfitCostFunction):
    """
    """

    def __init__(self, domain, data, maskZeroPotential = None,
                 conductivity_model=IPConductivityModelTemplate(),
                 surface_temperature = None, mask_surface_temperature = None,
                 mask_fixed_source = None,
                 sigma_src=None, pde_tol=1e-8, stationsFMT="e%s",
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, w1=1, reg_tol=None, logclip = 6,
                 thermal_conductivity=1, background_heat_source=None, surface_flux=None,
                 logger=None, **kargs):
        """
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.
        :param conductivity_model: conductivity model
        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
        :

        :param pde_tol: tolerance for solving the forward PDEs and recovering m from M.
        :param stationsFMT: format string to convert station id to mesh label
        :param m_ref: reference property function
        :param zero_mean_m: constrain m by zero mean.
        :param useLogMisfitDC: if set logarithm of DC data is used in misfit.
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolDC: relative tolerance for damping small DC data on misfit
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param length_scale: length scale, w0=(1/length scale)**2 is weighting factor for |grad m|^2 = |M|^2 in
                            the regularization term. If `None`, the term is dropped.
        :param theta: weigthing factor x-grad term
        :param fixTop: if set m[0]-m_ref[0] and m[1]-m_ref[1] are fixed at the top of the domain
                        rather than just on the left, right, front, back and bottom face.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param save_memory: if not set, the three PDEs for the three conponents of the inversion unknwon
                            with the different boundary conditions are held. Otherwise, boundary conditions
                            are updated for each component which requires refactorization which obviously
                            requires more time.
        :param logger: the logger, if not set, 'fingal.IPInversion.H2' is used.
        """
        if sigma_src == None:
            sigma_src = conductivity_model.sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversionBySource')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                         useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                         useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                         weightingMisfitDC=weightingMisfitDC,
                         logger=logger, **kargs)
        self.conductivity_model = conductivity_model
        self.logclip = logclip
        # PDE to recover temperature:
        if surface_temperature is None:
            surface_temperature = Scalar( self.conductivity_model.T_ref, Solution(domain))
        if mask_surface_temperature is None:
            z =Solution(domain).getX()[2]
            mask_surface_temperature = whereZero(z-sup(z))


        self.Tpde = setupERTPDE(self.domain)
        self.Tpde.getSolverOptions().setTolerance(pde_tol)
        #optionsG = self.Tpde.getSolverOptions()
        #from esys.escript.linearPDEs import SolverOptions
        #optionsG.setSolverMethod(SolverOptions.DIRECT)
        self.Tpde.setValue(A = thermal_conductivity * kronecker(3),
                           q=mask_surface_temperature, r=surface_temperature)
        ## get the initial temperature:
        if not surface_flux is None:
                self.Tpde.setValue( y=surface_flux)
        self.T_background = self.Tpde.getSolution()
        self.logger.debug("Background Temperature = %s" % (str(self.T_background)))
        self.Tpde.setValue(r=Data(), Y=Data(), y=Data())


        self.Q_background = background_heat_source
        self.Q_log = True
        # ... regularization
        if mask_fixed_source is None:
            x = domain.getX()
            qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
            qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])) )
            qz = whereZero(x[2] - inf(x[2]))
            self.mask_fixed_source  =   qx + qy + qz
        else:
            self.mask_fixed_source = wherePositive(mask_fixed_source)

        if not reg_tol:
            reg_tol = min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.Hpde.setValue(A=kronecker(3), q=self.mask_fixed_source)
        self.setW1(w1)
    #====
    def getTemperature(self, Q):
        """
        returns the temperature from property function M
        """
        self.Tpde.setValue( Y=Q, X=Data())
        T = self.Tpde.getSolution()
        self.logger.debug("Correction Temperature = %s" % (str(T)))
        T += self.T_background
        T = clip(T, minval=self.conductivity_model.T_min, maxval=self.conductivity_model.T_max)
        return T

    def setW1(self, w1):
        self.w1 = w1
        self.logger.debug(f'w1 = {self.w1:g}')

    def getSigma0(self, T, applyInterploation=False):
        """
        get the sigma_0 from temperature T
        """
        if applyInterploation:
            iT = interpolate(T, Function(self.domain))
        else:
            iT = T
        sigma_0 =  self.conductivity_model.getDCConductivity(T=iT)
        return sigma_0

    def getMn(self, T, applyInterploation=False):
        if applyInterploation:
            iT= interpolate(T, Function(self.domain))
        else:
            iT= T
        Mn  =  self.conductivity_model.getChargeability(T=iT)
        return Mn

    def getDsigma0DT(self, sigma_0, T):
        iT = interpolate(T, sigma_0.getFunctionSpace())
        return self.conductivity_model.getDsigma_0DT(T=iT)

    def getDMnDT(self, Mn, T):
        iT = interpolate(T, Mn.getFunctionSpace())
        return self.conductivity_model.getDMnDT(T=iT)

    def extractPropertyFunction(self, M):
        return M
    def getHeatSource(self, m, applyInterploation=False):
            """
            returns the heat source  from property function m
            """
            if hasattr(m, "getShape") and not m.getShape() == ():
                m = m[0]
            if applyInterploation:
                im = interpolate(m, Function(self.domain))
            else:
                im = m

            if self.Q_log:
                im = clip(im, minval=-self.logclip, maxval=self.logclip)
                Q = self.Q_background * exp(im)
            else:
                im = clip(im, minval=-exp(self.logclip), maxval=exp(self.logclip))
                Q = self.Q_background * (1+ im)
            return Q

    def getDHeatSourceDm(self, Q, m):
            """
            returns the derivative of heat source with respect of property function m
            """
            if self.Q_log:
                return Q
            else:
                return  self.Q_background


    def getArguments(self, M):
        """
        precalculation
        """
        m = self.extractPropertyFunction(M)
        im = interpolate(m, Function(self.domain))
        iQ = self.getHeatSource(im)
        T=self.getTemperature(iQ)
        iT = interpolate(T, Function(self.domain))
        iT_stations = self.grabValuesAtStations(T)
        self.logger.debug("Temperature = %s" % (str(iT)))
        self.logger.debug("Heat source = %s" % (str(iQ)))
        #self.logger.debug("Temperature at stations = %s" % (str(iT_stations)))

        isigma_0 = self.getSigma0(iT)
        isigma_0_stations = self.getSigma0(iT_stations)
        iMn = self.getMn(iT)
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return iQ, iT, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, M, iQ, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)
        gM = grad(M, where=iT.getFunctionSpace())
        R = 1. / 2. * integrate(self.w1 * length(gM) ** 2 )
        V = R + misfit_DC + misfit_IP
        self.logger.debug(
                f'misfit ERT, IP; reg, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e} = {V:e}')
        self.logger.debug(
                f'ratios ERT, IP; reg \t=  {misfit_DC / V * 100:g}, {misfit_IP / V * 100:g};  {R / V * 100:g}')
        return V

    #=======================================
    def getGradient(self, M, iQ, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gM = grad(M, where=iT.getFunctionSpace())
        X = self.w1 * gM

        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)
        Dsigma_0DT = self.getDsigma0DT(isigma_0, iT)
        DMnDT = self.getDMnDT(iMn, iT)
        Ystar = DMisfitDMn * DMnDT + DMisfitDsigma_0 * Dsigma_0DT
        self.Tpde.setValue(Y = Ystar, X=Data())
        Tstar = self.Tpde.getSolution()
        DmissfitDQ = interpolate(Tstar, where=X.getFunctionSpace())
        DQDm = self.getDHeatSourceDm(iQ, M)
        Y = DmissfitDQ * DQDm
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, M, iQ, iT, isigma_0, isigma_0_stations, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        p = self.Hpde.getSolution()/ self.w1
        self.logger.debug(f"search direction = {str(p)}.")
        return p

    def getDualProduct(self, M, r):
        """
        dual product of gradient `r` with increment `V`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0], M) + inner(r[1], grad(M)))

    def getNorm(self, M):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(M)

    def getSqueezeFactor(self, M, p):
        if self.Q_log:
            L = self.logclip
        else:
            L = exp(self.logclip)
        a=inf((self.logclip-sign(p)*M)/(abs(p)+1e-99) )
        if a > 0:
            return a
        else:
            return None

class InversionIPByTemperature(IPMisfitCostFunction):
    """
    """

    def __init__(self, domain, data, maskZeroPotential = None,
                 conductivity_model=IPConductivityModelTemplate(), logclip = 5,
                 surface_temperature = None, mask_surface_temperature = None,
                 sigma_src=None, pde_tol=1e-8, stationsFMT="e%s",
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, w1=1, reg_tol=None, temperature_variation_factor =0,
                 thermal_conductivity=1, heat_production=None, surface_flux=None,
                 logger=None, **kargs):
        """
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.
        :param conductivity_model: conductivity model
        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
        :

        :param pde_tol: tolerance for solving the forward PDEs and recovering m from M.
        :param stationsFMT: format string to convert station id to mesh label
        :param m_ref: reference property function
        :param zero_mean_m: constrain m by zero mean.
        :param useLogMisfitDC: if set logarithm of DC data is used in misfit.
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolDC: relative tolerance for damping small DC data on misfit
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param length_scale: length scale, w0=(1/length scale)**2 is weighting factor for |grad m|^2 = |M|^2 in
                            the regularization term. If `None`, the term is dropped.
        :param theta: weigthing factor x-grad term
        :param fixTop: if set m[0]-m_ref[0] and m[1]-m_ref[1] are fixed at the top of the domain
                        rather than just on the left, right, front, back and bottom face.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param save_memory: if not set, the three PDEs for the three conponents of the inversion unknwon
                            with the different boundary conditions are held. Otherwise, boundary conditions
                            are updated for each component which requires refactorization which obviously
                            requires more time.
        :param logger: the logger, if not set, 'fingal.IPInversion.H2' is used.
        """
        if sigma_src == None:
            sigma_src = conductivity_model.sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversionByTemperature')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                         useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                         useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                         weightingMisfitDC=weightingMisfitDC,
                         logger=logger, **kargs)
        self.logclip = logclip
        self.conductivity_model = conductivity_model
        # PDE to recover temperature:
        if surface_temperature is None:
            surface_temperature = Scalar( self.conductivity_model.T_ref, Solution(domain))
        if mask_surface_temperature is None:
            z =Solution(domain).getX()[2]
            mask_surface_temperature = whereZero(z-sup(z))

        self.Tpde = setupERTPDE(self.domain)
        self.Tpde.getSolverOptions().setTolerance(pde_tol)
        #optionsG = self.Tpde.getSolverOptions()
        #from esys.escript.linearPDEs import SolverOptions
        #optionsG.setSolverMethod(SolverOptions.DIRECT)
        self.Tpde.setValue(A = thermal_conductivity * kronecker(3),
                           q=mask_surface_temperature, r=surface_temperature)
        ## get the initial temperature:
        if not heat_production is None:
            self.Tpde.setValue( Y=heat_production)
        if not surface_flux is None:
                self.Tpde.setValue( y=surface_flux)
        self.T_background = self.Tpde.getSolution()
        self.logger.debug("Background Temperature = %s" % (str(self.T_background)))
        self.Tpde.setValue(r=Data(), Y=Data(), y=Data())

        x = domain.getX()
        qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
        qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])) )

        self.var_m_mask = wherePositive(interpolate(interpolate(mask_surface_temperature, ReducedFunctionOnBoundary(domain)), FunctionOnBoundary(domain)))
        self.var_m_b = temperature_variation_factor * self.var_m_mask

        if not reg_tol:
            reg_tol = min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.Hpde.setValue(A=kronecker(3) * w1, d = self.var_m_b , q=qx+qy)
        self.setW1(w1)
    #====
    def getTemperature(self, m, applyInterploation=False):
        """
        returns the temperature from property function M
        """
        if hasattr(m, "getShape") and not m.getShape() == ():
            m=m[0]
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        T = self.T_background * exp(im)
        return T

    def getDTemperatureDm(self, T, m):
        """
        returns the derivative of temperature with respect of property function m
        """
        return T

    def setW1(self, w1):
        self.w1 = w1
        self.Hpde.setValue(A=kronecker(3) * w1)
        self.logger.debug(f'w1 = {self.w1:g}')

    def getSigma0(self, T, applyInterploation=False):
        """
        get the sigma_0 from temperature T
        """
        if applyInterploation:
            iT = interpolate(T, Function(self.domain))
        else:
            iT = T
        sigma_0 =  self.conductivity_model.getDCConductivity(T=iT)
        return sigma_0

    def getMn(self, T, applyInterploation=False):
        if applyInterploation:
            iT= interpolate(T, Function(self.domain))
        else:
            iT= T
        Mn  =  self.conductivity_model.getChargeability(T=iT)
        return Mn

    def getDsigma0DT(self, sigma_0, T):
        iT = interpolate(T, sigma_0.getFunctionSpace())
        return self.conductivity_model.getDsigma_0DT(T=iT)

    def getDMnDT(self, Mn, T):
        iT = interpolate(T, Mn.getFunctionSpace())
        return self.conductivity_model.getDMnDT(T=iT)

    def extractPropertyFunction(self, m):
        return m

    def getArguments(self, m):
        """
        precalculation
        """
        #m = self.extractPropertyFunction(m)
        im = interpolate(m, Function(self.domain))
        iT=self.getTemperature(im)
        iT_stations = self.grabValuesAtStations(self.getTemperature(m))
        self.logger.debug("Temperature = %s" % (str(iT)))
        self.logger.debug("Temperature at stations = %s to %s " % (min(iT_stations), max(iT_stations)))

        isigma_0 = self.getSigma0(iT)
        isigma_0_stations = self.getSigma0(iT_stations)
        iMn = self.getMn(iT)
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return iT, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, m, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)
        gm = grad(m, where=iT.getFunctionSpace())
        R = 1. / 2. * integrate(self.w1 * length(gm) ** 2 )
        im_b=interpolate(m, FunctionOnBoundary(self.domain))
        Rbc= 1/2. *  integrate(self.var_m_b * im_b ** 2 )
        V = R + misfit_DC + misfit_IP + Rbc
        self.logger.debug(
                f'misfit ERT, IP; reg, vari, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e}; {Rbc:e}= {V:e}')
        self.logger.debug(
                f'ratios ERT, IP; reg \t=  {misfit_DC / V * 100:g}, {misfit_IP / V * 100:g};  {R / V * 100:g}')
        change_m = integrate(self.var_m_mask * im_b ** 2 )/integrate(self.var_m_mask)
        self.logger.debug(f"Change of property function on top boundary (mean) = {change_m:g}.")
        return V

    #=======================================
    def getGradient(self,m , iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gm = grad(m, where=iT.getFunctionSpace())

        X = self.w1 * gm

        im_b=interpolate(m, FunctionOnBoundary(self.domain))
        y = self.var_m_b * im_b

        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)
        Dsigma_0DT = self.getDsigma0DT(isigma_0, iT)
        DMnDT = self.getDMnDT(iMn, iT)
        DTdm = self.getDTemperatureDm(iT, m)
        Y = (DMisfitDMn * DMnDT + DMisfitDsigma_0 * Dsigma_0DT ) * DTdm
        return ArithmeticTuple(Y, X, y)

    def getInverseHessianApproximation(self, r, m, iT, isigma_0, isigma_0_stations, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0], y=r[2])
        p = self.Hpde.getSolution()
        self.logger.debug(f"search direction = {str(p)}.")
        return p

    def getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `V`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(r[0] * m + inner(r[1], grad(m))) + integrate( r[2] * m )

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)

    def getSqueezeFactor(self, m, p):
        a=inf((self.logclip-sign(p)*m)/(abs(p)+1e-99) )
        if a > 0:
            return a
        else:
            return None

class InversionIPByFluxWithWeakBC(IPMisfitCostFunction):
    """

    """

    def __init__(self, domain, data, maskZeroPotential = None,
                 conductivity_model=IPConductivityModelTemplate(),
                 set_temperature = 10., mask_set_temperature=Data(),
                 sigma_src=None, pde_tol=1e-8, stationsFMT="e%s", length_scale=None,
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, w1=1, reg_tol=None,
                 thermal_conductivity=1, heat_production=None, surface_flux=None,
                 mask_flux_variation = None, flux_variation_factor=0.,
                 logger=None, **kargs):
        """
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.
        :param conductivity_model: conductivity model
        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
        :

        :param pde_tol: tolerance for solving the forward PDEs and recovering m from M.
        :param stationsFMT: format string to convert station id to mesh label
        :param m_ref: reference property function
        :param zero_mean_m: constrain m by zero mean.
        :param useLogMisfitDC: if set logarithm of DC data is used in misfit.
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolDC: relative tolerance for damping small DC data on misfit
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param length_scale: length scale, w0=(1/length scale)**2 is weighting factor for |grad m|^2 = |M|^2 in
                            the regularization term. If `None`, the term is dropped.
        :param theta: weigthing factor x-grad term
        :param fixTop: if set m[0]-m_ref[0] and m[1]-m_ref[1] are fixed at the top of the domain
                        rather than just on the left, right, front, back and bottom face.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param save_memory: if not set, the three PDEs for the three conponents of the inversion unknwon
                            with the different boundary conditions are held. Otherwise, boundary conditions
                            are updated for each component which requires refactorization which obviously
                            requires more time.
        :param logger: the logger, if not set, 'fingal.IPInversion.H2' is used.
        """
        if sigma_src == None:
            sigma_src = conductivity_model.sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversionByTemperature')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                         useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                         useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                         weightingMisfitDC=weightingMisfitDC,
                         logger=logger, **kargs)
        self.conductivity_model = conductivity_model

        # PDE to set temperature:
        # if surface_temperature is None:
        #     surface_temperature = Scalar( self.conductivity_model.T_ref, Solution(domain))
        # if mask_surface_temperature is None:
        #     z =Solution(domain).getX()[2]
        #     mask_surface_temperature = whereZero(z-sup(z))
        #
        # x = domain.getX()
        # qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
        # qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])))
        # qz = whereZero(x[2] - inf(x[2]))

        self.Tpde = setupERTPDE(self.domain)
        self.Tpde.getSolverOptions().setTolerance(pde_tol)
        #optionsG = self.Tpde.getSolverOptions()
        #from esys.escript.linearPDEs import SolverOptions
        #optionsG.setSolverMethod(SolverOptions.DIRECT)
        self.Tpde.setValue(A = thermal_conductivity * kronecker(3),
                           q=mask_set_temperature, r=set_temperature)
        ## get the initial temperature:
        if not heat_production is None:
            self.Tpde.setValue( Y=heat_production )
        if not surface_flux is None:
                self.Tpde.setValue( y=surface_flux)
        self.T_background = self.Tpde.getSolution()
        self.logger.debug("Background Temperature = %s" % (str(self.T_background)))
        self.Tpde.setValue(r=Data(), Y=Data(), y=Data())
        saveSilo("T", T_bg= self.T_background, q=surface_flux, Tmask=mask_set_temperature)
        if length_scale is None:
            self.a = None
        else:
            self.a = (1. / length_scale) ** 2
        self.logger.info("Length scale factor is %s." % (str(length_scale)))
        # we allow variation of normal flux on all faces and the surface!
        self.surface = integrate(Scalar(1., FunctionOnBoundary(self.domain)))
        self.flux_variation_factor = flux_variation_factor
        self.scaled_flux_variation_factor = self.flux_variation_factor #/self.surface
        self.logger.info("Rescaled flux variation factor is %s." % (str(self.scaled_flux_variation_factor)))

        # ... regularization
        if not reg_tol:
            reg_tol = min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        # we assume that the flux correction terms are decoupling and use just the main diagonal of nn (which
        # is [0,0,1] for a flat surface). This to save compute time and memory when applying the preconditioner.
        self.Hpde = []
        nn = outer(self.domain.getNormal(), domain.getNormal())
        x = domain.getX()
        qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
        qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])))
        qz = whereZero(x[2] - inf(x[2]))

        #    # x = domain.getX()
        # qx = (whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0])))
        # qy = (whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1])))
        # qz = whereZero(x[2] - inf(x[2]))
        self.mask_flux_variation=mask_flux_variation
        if self.mask_flux_variation is None:
            qs = [qy + qz, qx + qz, qz]
            ff = 0.
        else:
            qs=[qy + qz, qx + qz, qx + qy ]
            ff = self.scaled_flux_variation_factor * self.mask_flux_variation

        for i, q in enumerate(qs):
            pde = setupERTPDE(self.domain)
            pde.getSolverOptions().setTolerance(reg_tol)
            pde.setValue(A=kronecker(3), d=nn[i, i] * ff, q = q)
            if self.a is not None:
                pde.setValue(D=self.a)
            self.Hpde.append(pde)
        self.setW1(w1)
    #====
    def getTemperature(self, M):
        """
        returns the temperature from property function M
        """
        self.Tpde.setValue(X=-M, Y =Data(), y=Data())
        T = self.Tpde.getSolution()
        self.logger.debug("Correction Temperature = %s" % (str(T)))
        T += self.T_background

        return T

    def setW1(self, w1):
        self.w1 = w1
        self.logger.debug(f'w1 = {self.w1:g}')

    def getSigma0(self, T, applyInterploation=False):
        """
        get the sigma_0 from temperature T
        """
        if applyInterploation:
            iT = interpolate(T, Function(self.domain))
        else:
            iT = T
        iT = clip(iT, minval=self.conductivity_model.T_min, maxval=self.conductivity_model.T_max)
        sigma_0 =  self.conductivity_model.getDCConductivity(T=iT)
        return sigma_0

    def getMn(self, T, applyInterploation=False):
        if applyInterploation:
            iT= interpolate(T, Function(self.domain))
        else:
            iT= T
        iT = clip(iT, minval=self.conductivity_model.T_min, maxval=self.conductivity_model.T_max)
        Mn  =  self.conductivity_model.getChargeability(T=iT)
        return Mn

    def getDsigma0DT(self, sigma_0, T):
        iT = interpolate(T, sigma_0.getFunctionSpace())
        Dsigma0DT = self.conductivity_model.getDsigma_0DT(T=iT)
        Dsigma0DT*=wherePositive(iT-self.conductivity_model.T_min) * wherePositive(self.conductivity_model.T_max-iT)
        return Dsigma0DT

    def getDMnDT(self, Mn, T):
        iT = interpolate(T, Mn.getFunctionSpace())
        DMnDT = self.conductivity_model.getDMnDT(T=iT)
        DMnDT *= wherePositive(iT - self.conductivity_model.T_min) * wherePositive(self.conductivity_model.T_max - iT)
        return DMnDT

    def extractPropertyFunction(self, M):
        return M

    def getArguments(self, M):
        """
        precalculation
        """
        M2 = self.extractPropertyFunction(M)
        iM = interpolate(M2, Function(self.domain))
        T=self.getTemperature(iM)
        iT = interpolate(T, Function(self.domain))
        iT_stations = self.grabValuesAtStations(T)
        self.logger.debug("Flux M[0] = %s" % (str(iM[0])))
        self.logger.debug("     M[1] = %s" % (str(iM[1])))
        self.logger.debug("     M[2] = %s" % (str(iM[2])))
        self.logger.debug("Temperature = %s" % (str(T)))
        #self.logger.debug("Temperature at stations = %s" % (str(iT_stations)))
        nnn=abs(inner(self.domain.getNormal(), M))
        S = integrate( self.mask_flux_variation * nnn**2)
        self.logger.debug(f'normal flux variation = {(S / self.surface)**0.5:g}, min = {inf(nnn+(1-self.mask_flux_variation)*1e99):g}, '
                          + f'max={sup(nnn - (1-self.mask_flux_variation) * 1e99):g}.')
        isigma_0 = self.getSigma0(iT)
        isigma_0_stations = self.getSigma0(iT_stations)
        iMn = self.getMn(iT)
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return iT, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, M, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        n=self.domain.getNormal()
        misfit_DC, misfit_IP = self.getMisfit(*args2)
        gM = grad(M, where=iT.getFunctionSpace())
        iM = interpolate(M, iT.getFunctionSpace())
        iiM = interpolate(M, n.getFunctionSpace())
        if self.a:
            R = 1. / 2. * integrate(self.w1 * (length(gM) ** 2 + self.a * length(iM) ** 2) )
        else:
            R = 1. / 2. * integrate(self.w1 * length(gM) ** 2  )

        S =  1. / 2. * integrate((self.w1 * self.scaled_flux_variation_factor) * inner(n, iiM) ** 2  )
        V = R + S + misfit_DC + misfit_IP
        self.logger.debug(
            f'       mis. ERT, mis. IP, reg., flux -> total \t=  {misfit_DC:g} + {misfit_IP:g} +  {R:g} + {S:g} = {V:g}')
        self.logger.debug(
            f'ratio: mis. ERT, mis. IP, reg., flux -> total \t=  {misfit_DC/V:e} + {misfit_IP/V:e} +  {R/V:e} + {S/V:e} = {1.:e}')

        return V

    #=======================================
    def getGradient(self, M, iT, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gM = grad(M, where=iT.getFunctionSpace())
        iM = interpolate(M, iT.getFunctionSpace())
        n = self.domain.getNormal()
        X = self.w1 * gM
        if self.a:
            Y = self.w1 * self.a * iM
        else:
            Y = Vector(0., iT.getFunctionSpace() )
        iiM = interpolate(M, n.getFunctionSpace())
        y = n * ( inner(n, iiM) * ( self.w1 * self.scaled_flux_variation_factor) )
        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)
        Dsigma_0DT = self.getDsigma0DT(isigma_0, iT)
        DMnDT = self.getDMnDT(iMn, iT)
        Ystar = DMisfitDMn * DMnDT + DMisfitDsigma_0 * Dsigma_0DT
        self.Tpde.setValue(Y = Ystar, X=Data(), y=Data())
        Tstar = self.Tpde.getSolution()
        Y += -grad(Tstar, where=Y.getFunctionSpace())
        return ArithmeticTuple(Y, X, y)

    def getInverseHessianApproximation(self, r, M, iT, isigma_0, isigma_0_stations, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        P = Data(0., (3,), Solution(self.domain))
        for k in [0, 1, 2]:
            self.Hpde[k].setValue(X=r[1][k], Y=r[0][k], y=r[2][k])
            P[k] = self.Hpde[k].getSolution() / self.w1
            self.logger.debug(f"search direction component {k} = {str(P[k])}.")
        return P

    def getDualProduct(self, M, r):
        """
        dual product of gradient `r` with increment `V`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        iiM=interpolate(M, r[2].getFunctionSpace())
        return integrate(inner(r[0], M) + inner(r[1], grad(M))) + integrate(inner(r[2], iiM))

    def getNorm(self, M):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(M)

    def getSqueezeFactor(self, M, p):
        return None
