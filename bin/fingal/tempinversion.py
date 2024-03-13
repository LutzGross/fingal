#
#
#
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
    def __init__(self, sigma_0_ref=1.51e-4, Mn_ref=5.05e-5, T_ref=15, P_0=2.9, P_Mn=2.5, **kargs):
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
