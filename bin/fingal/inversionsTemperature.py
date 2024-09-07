#
#
#
from .inversion import IPMisfitCostFunction
class IPConductivityModelTemplate(object):
    """
    this a template for class for providing an electric conductivity model for IP inversion
    the basic functionality is that it provides values for secondary conductivity sigma_0 as function of normalized chargeability M_n but may extend
    to other models such as temperature dependencies
    """

    def __init__(self, **kargs):
        pass

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
    def __init__(self, sigma_0_ref=1.52e-4, Mn_ref=5.054e-5, T_ref=15, P_0=2.9, P_Mn=2.5, **kargs):
        """
        :param T_ref: reference temperature
        :param sigma_0_ref: secondary conductivity at reference temperature
        :param Mn_ref: chargeability at at reference temperature
        :parem P_0: exponent of secondary conductivity  temperature model
        :parem P_Mn: exponent of chargeability temperature model
        """
        super().__init__(**kargs)
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

#=====================
class TemperatureBasedIPInversion(IPMisfitCostFunction):
    """
    Base class to run a IP inversion based on a temperature based model
    """
    def __init__(self, domain=None, data=None, weighting_misfit_ERT=0.5, pde_tol=1.e-8, stationsFMT="e%s",
                 conductivity_model=IPConductivityModelTemplate(), background_temperature=15.,
                 mask_fixed_temperature=None, useLogMisfitERT=False, useLogMisfitIP = False, logclip=5, **kargs):
        """
        :domain: pde domain
        :data: data, is `fingal.SurveyData`
        :sigma_0: reference conductivity
        :weighting_misfit_DC: weighting factor for
        :adjVs_ootStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        :param w0: weighting L2 regularization  int m^2
        :param w1: weighting H1 regularization  int grad(m)^2
        :mask_fixed_temperature_temperature: mask for fixed conductivity. needs to be set if w1>0 and w0=0
        """
        super().__init__(domain=domain, data=data, weighting_misfit_ERT=weighting_misfit_ERT,
                         sigma_background=sigma_background, useLogMisfitERT=useLogMisfitERT,
                         useLogMisfitERT=useLogMisfitIP,
                         pde_tol=pde_tol, stationsFMT=stationsFMT, logger=logger, **kargs)
        #==================TODO
        self.logclip = logclip
        self.EPS = EPSILON
        self.domain = domain
        self.logclip = logclip

        if mask_fixed_temperature is None:
            z = self.domain.getX()[2]
            self.mask_fixed_temperature = whereZero(z - sup(z))
        else:
            self.mask_fixed_temperature = mask_fixed_temperature
        self.conductivity_model = conductivity_model
        self.T_bg = background_temperature
        iTbg = interpolate(self.T_bg, where=Function(domain))
        sigma_0_bg = self.conductivity_model.getDCConductivity(T=self.T_bg)

        super().__init__(domain=domain, data=data, weighting_misfit_0=weighting_misfit_0, sigma_0_bg=sigma_0_bg,
                         pde_tol=pde_tol, stationsFMT=stationsFMT, **kargs)

        # regularization
        self.w2 = w2
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.setValue(A=self.w2 * kronecker(3), q=self.mask_fixed_temperature)
        self.Hpde.getSolverOptions().setTolerance(min(sqrt(pde_tol), 1e-3))

    def getTemperature(self, m, applyInterploation=False):
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        T = self.T_bg * exp(im)
        return T

    def getDTdm(self, T):
        return T

    def getArguments(self, m):
        """
        precalculated parameters:
        """
        # recover temperature (interpolated to elements)
        iT = self.getTemperature(m, applyInterploation=True)
        iMn = self.conductivity_model.getChargeability(T=iT)
        isigma_0 = self.conductivity_model.getDCConductivity(T=iT)
        txt2 = str(iT)
        txt1 = str(m)
        if getMPIRankWorld() == 0:
            self.logger.debug("m = " + txt1)
            self.logger.debug("temperature = " + txt2)

        args2 = self.getIPModelAndResponse(sigma_0=isigma_0, Mn=iMn)
        return iT, isigma_0, iMn, args2

    def getValue(self, m, *args):
        """
        return the value of the cost function
        """
        iT, isigma_0, iMn = args[0:3]

        misfit_2, misfit_0 = self.getMisfit(isigma_0, iMn, *args[3])

        H1 = self.w2 / 2 * integrate(length(grad(m)) ** 2)
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'reg: H1, misfits: ERT, IP =  {H1:e}, {misfit_0:e}, {misfit_2:e}')
        V = H1 + misfit_2 + misfit_0
        return V

    def getGradient(self, m, *args):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        iT, isigma_0, iMn = args[0:3]

        X = self.w2 * grad(m, where=iT.getFunctionSpace())
        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args[3])

        Dsigma_0DT = self.conductivity_model.getDsigma_0DT(iT)
        DMnDT = self.conductivity_model.getDMnDT(iT)
        DTDm = self.getDTdm(iT)
        Y = (DMisfitDMn * DMnDT + DMisfitDsigma_0 * Dsigma_0DT) * DTDm
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, m, *argsm, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        dm = self.Hpde.getSolution()

        txt = str(dm)
        if getMPIRankWorld() == 0:
            self.logger.debug(f"search direction = {txt}.")
        return dm

    def getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(r[0] * m + inner(r[1], grad(m)))

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)