"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException

from .tools import PoissonEquationZeroMean

from .tools import setupERTPDE, getSourcePotentials, makeMaskForOuterSurface, getAdditivePotentials, getSecondaryPotentials, DataMisfitQuad, DataMisfitLog
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple


class IP2MisfitCostFunction(CostFunction):
    """
    this the data misfit costfunction for IP

    Data are
       - the resistance as d_{ABMN}^{DC} = (U_A(X_M)-U_B(X_N)-U_A(X_M)+U_B(X_N))/I_{AB} for dipol injection AB of current I_{AB}
            and measurement over electrodes at X_M and X_N
       - over-voltage d_{ABMN}^IP = eta_{ABMN} * d_{ABMN}^{DC} with chargeability  eta_{ABMN}
    """

    def __init__(self, domain=None, data=None, sigma_0=1., pde_tol=1.e-8,
                 maskZeroPotential = None, stationsFMT="e%s", sigma_0_ref = None,
                 useLogMisfitIP=False, dataRTolIP=1e-4,
                 logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`. Resistance, charageability and errors are used.
        :param sigma_0: electric conductivity distribution.
        :type sigma_0: `Scalar` of type `Solution(domain)`.
        :param pde_tol: tolerance for the forward and adjoint PDE
        :param maskZeroPotential: mask of where potential is set to zero
        :param dataRTolIP: drop tolerance for small data, that is for data with values smaller than
                            dataRTolIP * max(resistance data).
        :param useLogMisfitIP: if True, logarithm of secondary potential data is used in the misfit
        :param dataRTolIP: drop tolerance for small data, that is for data with values smaller than
                            dataRToIP * max(secondary potential data).
        :stationsFMT: format string to map station keys A to mesh tags stationsFMT%A or None
        :param logger: `logging.Logger` if None, Logger `fingal` is used/created.
        """
        super().__init__(**kargs)
        self.domain = domain
        self.sigma_0=sigma_0
        if logger is None:
            self.logger = logging.getLogger('fingal.IP2Inversion')
        else:
            self.logger = logger
        if sigma_0_ref is None:
            self.sigma_0_ref = meanValue(self.sigma_0)
        else:
            self.sigma_0_ref = sigma_0_ref
        self.logger.debug("source & reference conductivity sigma_0_ref =  %s" % (str(self.sigma_0_ref)))
        if not maskZeroPotential or not Lsup(maskZeroPotential) > 0:
            raise ValueError("non-zero mask for zero potential must be given.")


        self.stationsFMT = stationsFMT
        self.maskZeroPotential = maskZeroPotential
        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.IP_pde = setupERTPDE(domain)
        self.IP_pde.getSolverOptions().setTolerance(pde_tol)
        self.IP_pde.setValue(q=self.maskZeroPotential)
        # self.DC_pde.getSolverOptions().setTrilinosParameter("reVs_ooe: type","full")
        self.data = data
        self.useLogMisfitIP=useLogMisfitIP
        self.dataRTolIP = dataRTolIP

        # extract values  by the Locator
        station_locations = []
        for k in data.getStationNumeration():
            station_locations.append(data.getStationLocationByKey(k))
        self.__grab_values_stations = Locator(DiracDeltaFunctions(domain), station_locations)
        self.setDCPotentials()
        # build the misfit data (indexed by source index):
        data_atol_IP= self.dataRTolIP * data.getMaximumSecondaryResistence()
        self.misfit_IP = {}  # secondary potentials
        nd_IP = 0 # counting number of data
        n_small_IP= 0 # number of dropped observations
        for A, B in self.data.injectionIterator():
            iA = self.data.getStationNumber(A)
            iB = self.data.getStationNumber(B)
            obs = self.data.getObservations(A, B)
            # ........... IP part .................................................................
            data_IP = np.array([self.data.getSecondaryResistenceData((A, B, M, N)) for M, N in obs])
            use_mask_IP =np.array( [ not self.data.isUndefined(v) for v in data_IP ] )
            data_IP = data_IP[use_mask_IP]
            obs_IP = [o for i, o in enumerate(obs) if use_mask_IP[i]]
            n_use_IP = len(data_IP)
            n_small_IP += np.count_nonzero(abs(data_IP) < data_atol_IP)
            if n_use_IP > 0:
                iMs = [self.data.getStationNumber(M) for M, N in obs_IP]
                iNs = [self.data.getStationNumber(N) for M, N in obs_IP]
                if self.useLogMisfitIP:
                    error_IP = np.array([self.data.getSecondaryResistenceRelError((A, B, M, N)) for M, N in obs_IP])  + self.dataRTolIP
                    self.misfit_IP[(iA, iB)] = DataMisfitLog(iMs=iMs, data=data_IP, iNs=iNs, injections=(A, B),
                                                             weightings=1. / error_IP ** 2 / n_use_IP)
                else:
                    #C = np.array([self.data.getChargeabilityData((A, B, M, N)) for M, N in obs_IP])
                    error_IP = np.array([self.data.getSecondaryResistenceError((A, B, M, N)) for M, N in obs_IP]) + data_atol_IP
                    self.misfit_IP[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_IP, iNs=iNs, injections=(A, B),
                                                              weightings=1. / error_IP ** 2 / n_use_IP)
                nd_IP += len(self.misfit_IP[(iA, iB)])
        self.logger.info(f"Data drop tolerance for secondary resistance is {data_atol_IP:e}.")
        self.logger.info(f"{nd_IP} IP data are records used. {n_small_IP} small values found.")
        if not nd_IP >0 :
            raise ValueError("No data for the inversion.")
        if nd_IP > 0:
            for iA, iB in self.misfit_IP:
                self.misfit_IP[(iA, iB)].rescaleWeight(1. / ( nd_IP) )
        self.ignoreIPMisfit(False)

    def ignoreIPMisfit(self, flag=False):
        """
        Switch on/off the IP misfit in the cost function. This is used for testing.
        """
        self.ignore_IPmisfit = flag
        if self.ignore_IPmisfit:
            self.logger.info(f"WARNING: **** IP Misfit is ignored in cost function ****")

    def grabValuesAtStations(self, m):
        """
        returns the values of m at the stations. The order is given by data.getStationNumeration()
        """
        return np.array(self.__grab_values_stations(m))

    def setDCPotentials(self):
        """
        return the primary electric potential for all injections (A,B) Vs_ooing conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_Field
        """
        sigma_src = self.sigma_0_ref
        source_potential = getSourcePotentials(self.domain, sigma_src, self.data, maskZeroPotential =self.maskZeroPotential, \
                                                    stationsFMT=self.stationsFMT, logger=self.logger)
        sigma_0_stations = self.grabValuesAtStations(self.sigma_0)
        additive_potentials_DC, src_potential_scale_DC = getAdditivePotentials(self.IP_pde,
                                                        sigma = self.sigma_0,
                                                        schedule = self.data,
                                                        sigma_stations = sigma_0_stations,
                                                        source_potential = source_potential,
                                                        sigma_src = sigma_src,
                                                        logger = self.logger)
        self.dc_potential_VA= {}
        self.dc_potential_VA_stations = {}
        for iA in additive_potentials_DC.keys():
            self.dc_potential_VA[iA] = src_potential_scale_DC[iA] * source_potential[iA] + additive_potentials_DC[iA]
            self.dc_potential_VA_stations[iA] = self.grabValuesAtStations(self.dc_potential_VA[iA])

        self.have_all_adjoints = set(self.data.getListofInjectionStationByIndex()).issubset(source_potential.keys())
        if self.have_all_adjoints:
            self.logger.info("Have all adjoint potentials.")

    def fitMnRef(self):
        """
        finds a new sigma_ref that gives a better data fit then
        """
        if len(self.misfit_IP) == 0:
            raise ValueError("No Data Available")
        beta = 1.
        a, b= 0.,0.
        if self.useLogMisfitIP:
            for iA, iB in self.misfit_IP:
                va = self.dc_potential_VA_stations[iA] - self.dc_potential_VA_stations[iB]
                dva = self.misfit_IP[(iA, iB)].getDifference(va)
                a+= sum(  self.misfit_IP[(iA, iB)].weightings * (self.misfit_IP[(iA, iB)].data - log(abs(dva)) ))
                b+= sum(  self.misfit_IP[(iA, iB)].weightings )
            if b>0:
                beta=exp(a/b)
        else:
            for iA, iB in self.misfit_IP:
                va = self.dc_potential_VA_stations[iA] - self.dc_potential_VA_stations[iB]
                dva = self.misfit_IP[(iA, iB)].getDifference(va)
                a+= sum(  self.misfit_IP[(iA, iB)].weightings * self.misfit_IP[(iA, iB)].data * dva )
                b+= sum(  self.misfit_IP[(iA, iB)].weightings * dva**2 )
            if b > 0:
                beta = a/b
        Mn_ref = beta / (1 - beta) * self.sigma_0
        return Mn_ref

    def getIPModelAndResponse(self, iMn):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and chargeaiBlity Mn
        they should be given on integration nodes.
        """
        isigma_0 = interpolate(self.sigma_0, iMn.getFunctionSpace())
        isigma_oo = isigma_0 + iMn
        self.logger.debug("Mn = %s" % (str(iMn)))
        self.logger.debug("sigma_oo = %s" % (str(isigma_oo)))
        secondary_potentials_IP = {}
        self.IP_pde.setValue(A=isigma_oo * kronecker(3), y_dirac=Data(), X=Data(), Y=Data())
        for iA in  self.dc_potential_VA:
            self.IP_pde.setValue(X=iMn * grad(self.dc_potential_VA[iA]))
            secondary_potentials_IP[iA] = self.IP_pde.getSolution()
        self.logger.info(f"{len(secondary_potentials_IP)} secondary IP potentials calculated.")
        secondary_potentials_IP_stations = {iA : self.grabValuesAtStations(secondary_potentials_IP[iA])
                                    for iA in secondary_potentials_IP}

        return secondary_potentials_IP, secondary_potentials_IP_stations

    def getMisfit(self, secondary_potentials_IP, secondary_potentials_IP_stations):
        """
        return the misfit in potential and chargeaiBlity weighted by weightingMisfitDC factor
        """
        misfit_IP = 0.
        for iA, iB in self.misfit_IP:
            misfit_IP += self.misfit_IP[(iA, iB)].getValue(
                secondary_potentials_IP_stations[iA] - secondary_potentials_IP_stations[iB])
        if self.ignore_IPmisfit:
            misfit_IP=0
        return misfit_IP

    def getDMisfit(self, iMn, secondary_potentials_IP, secondary_potentials_IP_stations):
        """
        returns the derivative of the misfit function with respect to sigma_0 and with respect to Mn
        """

        # ..... First we build the derivative of DC and IP misfit with respect to DC potential and secondary IP potenential
        SOURCES_IP = np.zeros((self.data.getNumStations(), self.data.getNumStations()), float)
        #potentials_DC_stations = {}
        #for iA in additive_potentials_DC_stations:
        #    potentials_DC_stations[iA] = src_potential_scale_DC[iA] * self.source_potentials_stations[iA] + additive_potentials_DC_stations[iA]

        DMisfitDMn = Scalar(0., iMn.getFunctionSpace() )
        for iA, iB in self.misfit_IP:
            dmis_IP = self.misfit_IP[(iA, iB)].getDerivative(
                secondary_potentials_IP_stations[iA] - secondary_potentials_IP_stations[iB])
            for i in range(len(self.misfit_IP[(iA, iB)])):
                iM = self.misfit_IP[(iA, iB)].iMs[i]
                iN = self.misfit_IP[(iA, iB)].iNs[i]
                SOURCES_IP[iA, iM] += dmis_IP[i]
                SOURCES_IP[iB, iM] -= dmis_IP[i]
                SOURCES_IP[iA, iN] -= dmis_IP[i]
                SOURCES_IP[iB, iN] += dmis_IP[i]

        if not self.ignore_IPmisfit:
            if self.have_all_adjoints:
                for iA in secondary_potentials_IP_stations.keys():
                    RA_star = Scalar(0., Solution(self.domain))
                    for M in self.data.getStationNumeration():
                        iM = self.data.getStationNumber(M)
                        VM = self.dc_potential_VA[iM]
                        WM = secondary_potentials_IP[iM]
                        RA_star+= SOURCES_IP[iA, iM] * (VM - WM)
                    # self.logger.debug("DC adjoint potential %d :%s" % (iA, str(VA_star)))
                    VA = self.dc_potential_VA[iA]
                    WA = secondary_potentials_IP[iA]
                    DMisfitDMn += inner(grad(RA_star), grad(VA - WA) )
            else:
                isigma_0 = interpolate(self.sigma_0, iMn.getFunctionSpace())
                isigma_oo = isigma_0 + iMn
                self.IP_pde.setValue(A=isigma_oo * kronecker(self.IP_pde.getDim()), X=Data(), Y=Data())
                DMisfitDMn = Scalar(0., self.DC_pde.getFunctionSpaceForCoefficient('Y'))
                for iA in secondary_potentials_IP_stations.keys():
                    VA = dc_potential_VA_stations[iA]
                    WA = secondary_potentials_IP[iA]
                    UA = VA - WA
                    s = Scalar(0., DiracDeltaFunctions(self.DC_pde.getDomain()))
                    for M in self.data.getStationNumeration():
                        iM = self.data.getStationNumber(M)
                        if self.stationsFMT is None:
                            s.setTaggedValue(M, SOURCES_IP[iA, iM])
                        else:
                            s.setTaggedValue(self.stationsFMT % M, SOURCES_IP[iA, iM])
                        self.IP_pde.setValue(y_dirac=s)
                    WA_star = self.IP_pde.getSolution()

                    DMisfitDMn += inner(grad(WA_star), grad(UA))
        # -------------------------
        self.logger.debug("%s adjoint potentials calculated." % len(secondary_potentials_IP_stations.keys()))
        return DMisfitDMn


class IP2InversionH1(IP2MisfitCostFunction):
    """
   Class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
   using H1-regularization on m-m_ref:

            m=log(Mn/Mn_ref)

   """

    def __init__(self, domain, data, maskZeroPotential, sigma_0 = 1., pde_tol= 1e-8,
                    stationsFMT="e%s", m_ref = None, fix_top = False, zero_mean_m = False,
                    useLogMisfitIP=False, dataRTolIP=1e-4, sigma_0_ref = None,
                    Mn_ref = 1.e-6,
                    w1=1, maskFixedProperty = None,
                    logclip=5, m_epsilon = 1e-18, reg_tol=None,
                    logger=None, **kargs):
        """       
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.

        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
        :param pde_tol: tolerance for solving the forward PDE
        :param stationsFMT: format string to convert station id to mesh label
        :param m_ref: reference property function
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity. If not set `meanValue` of `sigma_0` is used.
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param maskFixedProperty: mask of regions where m[0]-m_ref[0] and m[1]-m_ref[1] are fixed.
                                    if not set left, right, front, back and bottom are fixed.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param logger: the logger, if not set, 'fingal.IPInversion.H1' is used.
        """
        if logger is None:
            self.logger = logging.getLogger('fingal.IP2Inversion.H1')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, pde_tol= pde_tol, sigma_0 = sigma_0,
                            sigma_0_ref = sigma_0_ref,
                            maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                            useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                            logger=self.logger, **kargs)

        self.logclip = logclip
        self.m_epsilon = m_epsilon
        self.m_ref = m_ref
        self.zero_mean_m = zero_mean_m
        # its is assumed that sigma_ref and Mn_ref on the faces is not updated!!!
        if maskFixedProperty is None:
            x = self.IP_pde.getDomain().getX()
            qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
            qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
            qz = whereZero(x[2] - inf(x[2]))
            if fix_top:
                qz += whereZero(x[2] - sup(x[2]))
            self.mask_fixed_property  =   qx + qy + qz
        else:
            self.mask_fixed_property = wherePositive(maskFixedProperty)
        self.updateMnRef(Mn_ref)
        # regularization
        if reg_tol is None:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        if self.zero_mean_m:
            self.Hpde = PoissonEquationZeroMean(self.domain, A=kronecker(3)).setSecondaryCondition(Yh=1)
            self.Hpde.pde.getSolverOptions().setTolerance(reg_tol)
            self.logger.debug(f'Property function with zero mean.')
        else:
            self.Hpde = setupERTPDE(self.domain)
            self.Hpde.setValue(q=self.mask_fixed_property, A=kronecker(3))
            self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.setW1(w1)

    def setW1(self, w1):
        self.w1=w1
        self.logger.debug(f'w1 = {self.w1:g}')

    def updateMnRef(self, Mn_ref):
        """
        set a new reference conductivity
        """
        self.Mn_ref=Mn_ref
        self.logger.info("Reference normalised chargeability Mn_ref = %s" % (str(self.Mn_ref)))

    def getSigma0(self, m=None, applyInterploation=False):
        if applyInterploation:
            isigma_0 = interpolate(self.sigma_0, Function(self.domain))
        else:
            isigma_0 = self.sigma_0
        return isigma_0

    def getMn(self, m, applyInterploation=False):
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        Mn = self.Mn_ref * exp(im)
        return Mn

    def getDMnDm(self, Mn, m):
        return Mn

    def extractPropertyFunction(self, dm):
        if self.m_ref:
            m = dm + self.m_ref
        else:
            m = dm
        return m

    def getArguments(self, dm):
        """
        precalculation
        """
        m = self.extractPropertyFunction(dm)
        im = interpolate(m, Function(self.domain))
        #im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % (str(im)))
        iMn = self.getMn(im)
        
        args2 = self.getIPModelAndResponse(iMn)
        return im, iMn, args2

    def getValue(self, dm, im, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_IP = self.getMisfit(*args2)

        gdm=grad(dm, where=im.getFunctionSpace())
        R = self.w1 / 2 * integrate(length(gdm) ** 2 )

        V =  R + misfit_IP
        self.logger.debug(
                f'misfit IP, reg = total \t;  {misfit_IP:e};  {R:e} = {V:e}')
        self.logger.debug(
                f'ratios IP, reg [%] \t = {misfit_IP/V*100:g};  {R/V*100:g}')
        return V

    def getGradient(self, dm,  im, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gdm=grad(dm, where=im.getFunctionSpace())
        X = self.w1 * gdm
        DMisfitDMn  = self.getDMisfit(iMn, *args2)
        DMnDm = self.getDMnDm(iMn, im)
        Y = DMisfitDMn * DMnDm

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, dm, iMn, *args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        if self.zero_mean_m:
                dm = self.Hpde.getSolution(X=r[1], Y=r[0])/self.w1
        else:
                self.Hpde.setValue(X=r[1], Y=r[0])
                dm = self.Hpde.getSolution()/self.w1
        self.logger.debug(f"search direction = {str(dm)}.")
        return dm

    def getDualProduct(self, dm, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0], dm) + inner(r[1], grad(dm)))

    def getNorm(self, dm):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(dm)

    def getSqueezeFactor(self, dm, p):
        if self.m_ref:
            m = dm + self.m_ref
        else:
            m = dm

        a=inf((self.logclip-sign(p)*m)/(abs(p)+1e-99) )
        if a > 0:
            return a
        else:
            return None


class IP2InversionH2(IP2MisfitCostFunction):
    """
   Class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
   using H2-regularization on M = grad(m-m_ref) (as 3 components):

        m=log(Mn/Mn_ref)

   """

    def __init__(self, domain, data, maskZeroPotential, sigma_src=None, pde_tol=1e-8,
                 stationsFMT="e%s", m_ref=None, length_scale=None,
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, sigma_0_ref=1e-4, Mn_ref=1.e-6,
                 w1=1, fixTop=False,  zero_mean_m = True,
                 logclip=5, m_epsilon=1e-18, reg_tol=None, save_memory=False,
                 logger=None, **kargs):
        """
        :param domain: PDE & inversion domain
        :param data: survey data, is `fingal.SurveyData`. Resistence and secondary potential data are required.
        :param maskZeroPotential: mask of locations where electric potential is set to zero.
        :param sigma_src: background conductivity, used to calculate the source potentials. If not set sigma_0_ref
        is used.
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
            sigma_src = sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversion.H2')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                         useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                         useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                         weightingMisfitDC=weightingMisfitDC,
                         logger=logger, **kargs)

        self.logclip = logclip
        self.m_epsilon = m_epsilon
        self.m_ref = m_ref
        self.zero_mean_m =  zero_mean_m
        # its is assumed that sigma_ref and Mn_ref on the faces is not updated!!!
        x=self.domain.getX()
        qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
        qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
        qz = whereZero(x[2] - inf(x[2]))
        if fixTop:
            qz += whereZero(x[2] - sup(x[2]))

        self.updateSigma0Ref(sigma_0_ref)
        self.updateMnRef(Mn_ref)
        # regularization
        if self.zero_mean_m:
            self.logger.info(f'Property function with zero mean.')
        else:
            if fixTop:
                self.logger.info(f'Property function zero on all boundaries.')
            else:
                self.logger.info(f'Property function free on top and zero on all other boundaries.')

        if length_scale is None:
            self.a = 0.
        else:
            self.a = ( 1./length_scale )**2
        self.logger.info("Length scale factor is %s." % (str(length_scale)))

        # recovery of propery function from M:
        if self.zero_mean_m:
            self.mpde = PoissonEquationZeroMean(self.domain).setSecondaryCondition(Yh=1)
            self.logger.debug(f'Property function with zero mean.')
        else:
            self.mpde = setupERTPDE(self.domain)
            self.mpde.getSolverOptions().setTolerance(pde_tol)
            self.mpde.setValue(A=kronecker(3), q= qx + qy + qz )

        # regularization
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)

        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.save_memory = save_memory
        if self.save_memory:
            self.Hpde = setupERTPDE(self.domain)
            self.Hpde.getSolverOptions().setTolerance(reg_tol)
            self.Hpde_qs = [ qy + qz, qx + qz,  qx + qy ]
            self.Hpde.setValue(A=kronecker(3), D=self.a)
        else:
            self.Hpdes = []
            for k, q in enumerate( [ qy + qz, qx + qz,  qx + qy ]):
                pde =  setupERTPDE(self.domain)
                pde.getSolverOptions().setTolerance(reg_tol)
                pde.setValue(q=q, A=kronecker(3), D=self.a)
                self.Hpdes.append(pde)


        self.setW1(w1)

    def setW1(self, w1):
        self.w1 = w1
        self.logger.debug(f'w1 = {self.w1:g}')

    def updateMnRef(self, Mn_ref):
        """
        set a new reference conductivity
        """
        self.Mn_ref = Mn_ref
        self.logger.info("Reference normalised chargeability Mn_ref = %s" % (str(self.Mn_ref)))

    def getMn(self, m, applyInterploation=False):
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        Mn = self.Mn_ref * exp(im)
        return Mn

    def getDMnDm(self, Mn, m):
        return Mn

    def extractPropertyFunction(self, dM):
        if self.zero_mean_m:
            dm = self.mpde.getSolution(X=interpolate(dM, Function(dM.getDomain())))
        else:
            self.mpde.setValue(X=interpolate(dM, Function(dM.getDomain())), Y=Data(), y=Data())
            dm=self.mpde.getSolution()
        if self.m_ref:
            m=dm+self.m_ref
        else:
            m=dm
        return m

    def getArguments(self, dM):
        """
        precalculation
        """
        m = self.extractPropertyFunction(dM)
        im = interpolate(m, Function(self.domain))
        self.logger.debug("m = %s" % (str(im)))
        iMn = self.getMn(im)

        args2 = self.getIPModelAndResponse(iMn)
        return im, iMn, args2

    def getValue(self, dM, im, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_IP = self.getMisfit(*args2)

        gdM = grad(dM, where=im.getFunctionSpace())
        idM = interpolate(dM, im.getFunctionSpace() )
        R = self.w1 / 2. * integrate( length(gdM) ** 2 + self.a * length(idM)**2 )
        V = R + misfit_IP
        self.logger.debug(
                f'misfit IP; reg, total \t=  {misfit_IP:e};  {R:e} = {V:e}')
        self.logger.debug(
                f'ratios IP; reg \t=  {misfit_IP / V * 100:g};  {R / V * 100:g}')

        return V

    def getGradient(self, dM, im, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gdM = grad(dM, where=im.getFunctionSpace())
        idM = interpolate(dM, im.getFunctionSpace() )
        X = self.w1 * gdM
        Y = self.w1 * self.a * idM

        DMisfitDMn = self.getDMisfit(iMn, *args2)
        DMnDm = self.getDMnDm(iMn, im)
        Ystar= DMisfitDMn * DMnDm
        if self.zero_mean_m:
            dmstar = self.mpde.getSolution( Y = Ystar )
        else:
            self.mpde.setValue(X=Data(), Y = Ystar)
            dmstar=self.mpde.getSolution()
        Y+= grad(dmstar, where=X.getFunctionSpace())

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, dM, isigma_0, isigma_0_stations,
                                       iMn, *args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """

        P = Data(0., (6,), Solution(self.domain))
        for k in [0, 1, 2]:
            if self.save_memory:
                self.Hpde.setValue(q=self.Hpde_qs[k])
            for i in [0, 1]:
                if self.save_memory:
                    self.Hpde.setValue(X=r[1][i*3+k], Y=r[0][i*3+k])
                    P[i * 3 + k] = self.Hpde.getSolution() / self.w1
                else:
                    self.Hpdes[k].setValue(X=r[1][i * 3 + k], Y=r[0][i * 3 + k])
                    P[i*3+k] = self.Hpdes[k].getSolution() / self.w1[i]
                self.logger.debug(f"search direction component {i*3+k} = {str(P[i*3+k])}.")
        return P

    def getDualProduct(self, dM, r):
        """
        dual product of gradient `r` with increment `V`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0], dM) + inner(r[1], grad(dM)))

    def getNorm(self, dM):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(dM)

    def getSqueezeFactor(self, dM, p):
            return None
