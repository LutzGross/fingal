"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException


from .tools import setupERTPDE, getSourcePotentials, makeMaskForOuterSurface, getAdditivePotentials, getSecondaryPotentials, DataMisfitQuad, DataMisfitLog
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple


class IPMisfitCostFunction(CostFunction):
    """
    this the data misfit costfunction for IP

    Data are
       - the resistance as d_{ABMN}^{DC} = (U_A(X_M)-U_B(X_N)-U_A(X_M)+U_B(X_N))/I_{AB} for dipol injection AB of current I_{AB}
            and measurement over electrodes at X_M and X_N
       - over-voltage d_{ABMN}^IP = eta_{ABMN} * d_{ABMN}^{DC} with chargeability  eta_{ABMN}
    """

    def __init__(self, domain=None, data=None, sigma_src=1., pde_tol=1.e-8,
                 maskZeroPotential = None, stationsFMT="e%s",
                 useLogMisfitDC=False, dataRTolDC=1e-4,
                 useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC = 1,
                 logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`. Resistance, charageability and errors are used.
        :param sigma_src: electric conductivity used to calculate the source potential.
        :param pde_tol: tolerance for the forward and adjoint PDE
        :param maskZeroPotential: mask of where potential is set to zero
        :param useLogMisfitDC: if True, logarithm of data is used in the misfit
        :param dataRTolDC: drop tolerance for small data, that is for data with values smaller than
                            dataRTolDC * max(resistance data).
        :param useLogMisfitIP: if True, logarithm of secondary potential data is used in the misfit
        :param dataRTolIP: drop tolerance for small data, that is for data with values smaller than
                            dataRToIP * max(secondary potential data).
        :stationsFMT: format string to map station keys A to mesh tags stationsFMT%A or None
        :param logger: `logging.Logger` if None, Logger `fingal` is used/created.
        """
        super().__init__(**kargs)
        if logger is None:
            self.logger = logging.getLogger('fingal')
        else:
            self.logger = logger
        if not maskZeroPotential or not Lsup(maskZeroPotential) > 0:
            raise ValueError("non-zero mask for zero potential must be given.")

        self.domain = domain
        self.stationsFMT = stationsFMT
        self.maskZeroPotential = maskZeroPotential
        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.DC_pde = setupERTPDE(domain)
        self.DC_pde.getSolverOptions().setTolerance(pde_tol)
        self.DC_pde.setValue(q=self.maskZeroPotential)
        self.IP_pde = setupERTPDE(domain)
        self.IP_pde.getSolverOptions().setTolerance(pde_tol)
        self.IP_pde.setValue(q=self.maskZeroPotential)
        # self.DC_pde.getSolverOptions().setTrilinosParameter("reVs_ooe: type","full")
        self.data = data
        self.useLogMisfitDC=useLogMisfitDC
        self.dataRTolDC = dataRTolDC
        self.useLogMisfitIP=useLogMisfitIP
        self.dataRTolIP = dataRTolIP
        self.setMisfitWeighting(weightingMisfitDC)

        # extract values  by the Locator
        station_locations = []
        for k in data.getStationNumeration():
            station_locations.append(data.getStationLocationByKey(k))
        self.__grab_values_stations = Locator(DiracDeltaFunctions(domain), station_locations)
        # self.grabValues=Locator(Solution(domain), station_locations)
        self.sigma_src = sigma_src  # needs to be a constant.
        self.setSourcePotentials()  # S_s is in the notes
        # build the misfit data (indexed by source index):
        data_atol_DC= self.dataRTolDC * data.getMaximumResistence()
        data_atol_IP= self.dataRTolIP * data.getMaximumSecondaryResistence()

        self.misfit_IP = {}  # secondary potentials
        self.misfit_DC = {}  #
        nd_DC, nd_IP = 0,0 # counting number of data
        n_small_DC, n_small_IP= 0,0 # number of dropped observations
        for A, B in self.data.injectionIterator():
            iA = self.data.getStationNumber(A)
            iB = self.data.getStationNumber(B)
            obs = self.data.getObservations(A, B)
            # ........... DC part .................................................................
            data_DC = np.array([self.data.getResistenceData((A, B, M, N)) for M, N in obs])
            obs_DC = obs
            n_small_DC += np.count_nonzero(abs(data_DC) < data_atol_DC)
            n_use_DC=len(data_DC)
            if n_use_DC > 0:
                iMs = [self.data.getStationNumber(M) for M, N in obs_DC]
                iNs = [self.data.getStationNumber(N) for M, N in obs_DC]
                if self.useLogMisfitDC:
                    error_DC = np.array([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs_DC]) + self.dataRTolDC
                    self.misfit_DC[(iA, iB)] =  DataMisfitLog(iMs=iMs, data=data_DC, iNs=iNs, injections=(A,B), weightings=1. / error_DC**2/n_use_DC)
                else:
                    error_DC = np.array([self.data.getResistenceError((A, B, M, N)) for M, N in obs_DC]) + data_atol_DC
                    self.misfit_DC[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B), weightings =1./error_DC**2/n_use_DC)
                nd_DC+= n_use_DC
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
                    error_IP = np.array([self.data.getSecondaryResistenceError((A, B, M, N)) for M, N in obs_IP]) + data_atol_IP
                    self.misfit_IP[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_IP, iNs=iNs, injections=(A, B),
                                                              weightings=1. / error_IP ** 2 / n_use_IP)
                nd_IP += len(self.misfit_IP[(iA, iB)])
        self.logger.info(f"Data drop tolerance for resistance is {data_atol_DC:e}.")
        self.logger.info(f"{nd_DC} DC data are records used. {n_small_DC} small values found.")
        self.logger.info(f"Data drop tolerance for secondary resistance is {data_atol_IP:e}.")
        self.logger.info(f"{nd_IP} IP data are records used. {n_small_IP} small values found.")
        if not nd_IP + nd_DC >0 :
            raise ValueError("No data for the inversion.")
        if nd_DC > 0:
            for iA, iB in self.misfit_DC:
                self.misfit_DC[(iA, iB)].rescaleWeight(1. / (2 * nd_DC) )
        if nd_IP > 0:
            for iA, iB in self.misfit_IP:
                self.misfit_IP[(iA, iB)].rescaleWeight(1. / (2 * nd_IP))
        self.ignoreERTMisfit(False)
        self.ignoreIPMisfit(False)

    def ignoreERTMisfit(self, flag=False):
        """
        Switch on/off the ERT misfit in the cost function. This is used for testing.
        """
        self.ignore_ERTmisfit = flag
        if self.ignore_ERTmisfit:
            self.logger.info(f"WARNING: **** ERT Misfit is ignored in cost function ****")
    def ignoreIPMisfit(self, flag=False):
        """
        Switch on/off the IP misfit in the cost function. This is used for testing.
        """
        self.ignore_IPmisfit = flag
        if self.ignore_IPmisfit:
            self.logger.info(f"WARNING: **** IP Misfit is ignored in cost function ****")

    def setMisfitWeighting(self, weightingMisfitDC=1):
        """
        sets the relative weighting of DC and IP part in the misfit.
        """
        assert weightingMisfitDC>=0
        self.weightingMisfit_DC = weightingMisfitDC/(weightingMisfitDC+1)
        self.weightingMisfit_IP = max(0., 1 - self.weightingMisfit_DC)
        self.logger.info(f"Weighting of DC & IP in misfit = {self.weightingMisfit_DC} + {self.weightingMisfit_IP}")

    def grabValuesAtStations(self, m):
        """
        returns the values of m at the stations. The order is given by data.getStationNumeration()
        """
        return np.array(self.__grab_values_stations(m))

    def setSourcePotentials(self):
        """
        return the primary electric potential for all injections (A,B) Vs_ooing conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_Field
        """
        self.logger.debug("source conductivity sigma_src =  %s" % (str(self.sigma_src)))
        self.source_potential = getSourcePotentials(self.domain, self.sigma_src, self.data, maskZeroPotential =self.maskZeroPotential, \
                                                    stationsFMT=self.stationsFMT, logger=self.logger)
        self.source_potentials_stations = { iA : self.grabValuesAtStations(self.source_potential[iA])
                                               for iA in self.source_potential }

    def getIPModelAndResponse(self,sigma_0, sigma_0_stations, Mn):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and chargeaiBlity Mn
        they should be given on integration nodes.
        """
        sigma_oo = sigma_0 + Mn

        self.logger.debug("sigma_0 = %s" % (str(sigma_0)))
        self.logger.debug("Mn = %s" % (str(Mn)))
        self.logger.debug("sigma_oo = %s" % (str(sigma_oo)))
        self.logger.debug("sigma_0_stations = %g - %g " % (min(sigma_0_stations), max(sigma_0_stations)))

        additive_potentials_DC, src_potential_scale_DC = getAdditivePotentials(self.DC_pde,
                                                        sigma = sigma_0,
                                                        schedule = self.data,
                                                        sigma_stations = sigma_0_stations,
                                                        source_potential = self.source_potential,
                                                        sigma_src = self.sigma_src,
                                                        logger = self.logger)
        additive_potentials_DC_stations = {iA : self.grabValuesAtStations(additive_potentials_DC[iA])
                                               for iA in additive_potentials_DC}


        secondary_potentials_IP = getSecondaryPotentials(self.IP_pde, sigma_oo,
                                                Mn, source_potential=self.source_potential,
                                                additive_potential_DC=additive_potentials_DC,
                                                src_potential_scale_DC= src_potential_scale_DC,
                                                logger=self.logger)
        secondary_potentials_IP_stations = {iA : self.grabValuesAtStations(secondary_potentials_IP[iA])
                                    for iA in secondary_potentials_IP}

        return additive_potentials_DC, src_potential_scale_DC, additive_potentials_DC_stations, secondary_potentials_IP, secondary_potentials_IP_stations

    def getMisfit(self, additive_potentials_DC, src_potential_scale_DC, additive_potentials_DC_stations, secondary_potentials_IP, secondary_potentials_IP_stations, *args):
        """
        return the misfit in potential and chargeaiBlity weighted by weightingMisfitDC factor
        """
        potentials_DC_stations = {}
        for iA in additive_potentials_DC_stations:
            potentials_DC_stations[iA] = src_potential_scale_DC[iA] * self.source_potentials_stations[iA] + additive_potentials_DC_stations[iA]

        misfit_IP, misfit_DC = 0., 0.
        for iA, iB in self.misfit_IP:
            misfit_IP += self.misfit_IP[(iA, iB)].getValue(
                secondary_potentials_IP_stations[iA] - secondary_potentials_IP_stations[iB])
            misfit_DC += self.misfit_DC[(iA, iB)].getValue(
                 potentials_DC_stations[iA] - potentials_DC_stations[iB])

        if self.ignore_ERTmisfit:
            misfit_DC=0
        if self.ignore_IPmisfit:
            misfit_IP=0

        return misfit_DC * self.weightingMisfit_DC, misfit_IP * self.weightingMisfit_IP

    def getDMisfit(self, sigma_0, Mn, additive_potentials_DC, src_potential_scale_DC, additive_potentials_DC_stations, secondary_potentials_IP, secondary_potentials_IP_stations, *args):
        """
        returns the derivative of the misfit function with respect to sigma_0 and with respect to Mn
        """

        # ..... First we build the derivative of DC and IP misfit with respect to DC potential and secondary IP potenential
        SOURCES_DC = np.zeros( (self.data.getNumStations(), self.data.getNumStations()), float)
        SOURCES_IP = np.zeros((self.data.getNumStations(), self.data.getNumStations()), float)
        potentials_DC_stations = {}
        for iA in additive_potentials_DC_stations:
            potentials_DC_stations[iA] = src_potential_scale_DC[iA] * self.source_potentials_stations[iA] + additive_potentials_DC_stations[iA]
        if not self.ignore_ERTmisfit:
            for iA, iB in self.misfit_DC:
                dmis_DC = self.misfit_DC[(iA, iB)].getDerivative(
                    potentials_DC_stations[iA] - potentials_DC_stations[iB])
                for i in range(len(self.misfit_DC[(iA, iB)])):
                    iM = self.misfit_DC[(iA, iB)].iMs[i]
                    iN = self.misfit_DC[(iA, iB)].iNs[i]
                    #ABMN
                    SOURCES_DC[iA, iM] += dmis_DC[i]
                    SOURCES_DC[iB, iM] -= dmis_DC[i]
                    SOURCES_DC[iA, iN] -= dmis_DC[i]
                    SOURCES_DC[iB, iN] += dmis_DC[i]
            SOURCES_DC*=self.weightingMisfit_DC

        if not self.ignore_IPmisfit:
            for iA, iB in self.misfit_IP:
                dmis_IP = self.misfit_IP[(iA, iB)].getDerivative(
                    secondary_potentials_IP_stations[iA] - secondary_potentials_IP_stations[iB])
                for i in range(len(self.misfit_IP[(iA, iB)])):
                    iM = self.misfit_IP[(iA, iB)].iMs[i]
                    iN = self.misfit_IP[(iA, iB)].iNs[i]
                    #ABMN
                    SOURCES_IP[iA, iM] += dmis_IP[i]
                    SOURCES_IP[iB, iM] -= dmis_IP[i]
                    SOURCES_IP[iA, iN] -= dmis_IP[i]
                    SOURCES_IP[iB, iN] += dmis_IP[i]
            SOURCES_IP *= self.weightingMisfit_IP
        sigma_oo = (sigma_0 + Mn)
        self.DC_pde.setValue(A=sigma_0 * kronecker(self.DC_pde.getDim()), X=Data(), Y=Data())
        self.IP_pde.setValue(A=sigma_oo * kronecker(self.IP_pde.getDim()), X=Data(), Y=Data())
        DMisfitDsigma_0 = Scalar(0., self.DC_pde.getFunctionSpaceForCoefficient('Y'))
        DMisfitDMn = Scalar(0., self.DC_pde.getFunctionSpaceForCoefficient('Y'))
        for iA in additive_potentials_DC.keys():
            VA = src_potential_scale_DC[iA] * self.source_potential[iA] + additive_potentials_DC[iA]
            WA = secondary_potentials_IP[iA]
            s = Scalar(0., DiracDeltaFunctions(self.DC_pde.getDomain()))
            for M in self.data.getStationNumeration():
                iM = self.data.getStationNumber(M)
                if self.stationsFMT is None:
                    s.setTaggedValue(M, SOURCES_IP[iA, iM])
                else:
                    s.setTaggedValue(self.stationsFMT % M, SOURCES_IP[iA, iM])
            self.IP_pde.setValue(y_dirac=s)
            WA_star = self.IP_pde.getSolution()
            s = Scalar(0., DiracDeltaFunctions(self.DC_pde.getDomain()))
            for M in self.data.getStationNumeration():
                iM = self.data.getStationNumber(M)
                if self.stationsFMT is None:
                    s.setTaggedValue(M, SOURCES_DC[iA, iM])
                else:
                    s.setTaggedValue(self.stationsFMT % M, SOURCES_DC[iA, iM])
            self.DC_pde.setValue(y_dirac=s, X=-Mn * grad(WA_star))
            VA_star = self.DC_pde.getSolution()
            UA = WA + VA
            DMisfitDsigma_0 -= inner(grad(VA_star), grad(VA)) + inner(grad(WA_star), grad(WA))
            DMisfitDMn -= inner(grad(WA_star), grad(UA))
        # -------------------------
        self.logger.debug("%s adjoint potentials calculated." % len(additive_potentials_DC.keys()))
        return DMisfitDsigma_0,  DMisfitDMn


class IPInversionH1(IPMisfitCostFunction):
    """
    Base class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
    """

    def __init__(self,  domain=None, data=None, sigma_src=None, pde_tol= 1e-8,
                    maskZeroPotential=None, stationsFMT="e%s",
                    useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                    weightingMisfitDC=1, sigma_0_ref=1e-4, Mn_ref = 1.e-6,
                    w1=1, theta=0., maskFixedProperty = None,
                    logclip=5, m_epsilon = 1e-18, reg_tol=None,
                    logger=None, **kargs):
        """       
        :domain: pde domain
        :data: survey data, is `fingal.SurveyData`
        :maskFixedProperty: mask of region where Mn and sigma_0 are fixed (set to reference values).
            If not set the bottom of the domain is used.
        :sigma_background: background conductivity, used to calculate the injection potentials
        :sigma_0_ref: reference DC conductivity
        :Mn_ref: reference normalized chargeability (must be positive)
        :param w1: weighting H1 regularization  int grad(m)^2
        :param theta: weighting cross-gradient term
        :weightingMisfitDC: weighting factor for ERT part cost function, use 0.5
        :adjVs_ootStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format Vs_ooed to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        :EPSILON: absolute tolerance for zero values.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol= pde_tol,
                            maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                            useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                            useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                            weightingMisfitDC=weightingMisfitDC,
                            logger=logger, **kargs)

        self.logclip = logclip
        self.m_epsilon = m_epsilon
        # its is assumed that sigma_ref and Mn_ref on the faces is not updated!!!
        if maskFixedProperty is None:
            x = self.DC_pde.getDomain().getX()
            qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
            qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
            qz = whereZero(x[2] - inf(x[2]))
            self.mask_fixed_property  =   qx + qy + qz
        else:
            self.mask_fixed_property = wherePositive(maskFixedProperty)
        self.updateSigma0Ref(sigma_0_ref)
        self.updateMnRef(Mn_ref)
        # regularization

        if reg_tol is None:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.setValue(q=self.mask_fixed_property)
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.setW1andTheta(w1, theta)



    def setW1andTheta(self, w1, theta):
        self.w1=w1
        self.theta = theta
        self.Hpde.setValue(A=(self.w1+self.theta) * kronecker(3))
        self.logger.debug(f'w1 = {self.w1:g}')
        self.logger.debug(f'theta = {self.theta:g}')


    def updateSigma0Ref(self, sigma_0_ref):
        """
        set a new reference conductivity
        """
        self.sigma_0_ref=sigma_0_ref
        self.logger.info("Reference conductivity sigma_0_ref is set to %s" % (str(self.sigma_0_ref)))
    def updateMnRef(self, Mn_ref):
        """
        set a new reference conductivity
        """
        self.Mn_ref=Mn_ref
        self.logger.info("Reference normalised chargeability Mn_ref = %s" % (str(self.Mn_ref)))

    def getSigma0(self, m, applyInterploation=False):
        if hasattr(m, "getShape") and not m.getShape() == ():
            m=m[0]
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        sigma_0 = self.sigma_0_ref * exp(im)
        return sigma_0

    def getDsigma0Dm(self, sigma_0, m):
        return sigma_0

    def getMn(self, m, applyInterploation=False):
        if hasattr(m, "getShape") and not m.getShape() == ():
            m=m[1]
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        Mn = self.Mn_ref * exp(im)
        return Mn

    def getDMnDm(self, Mn, m):
        return Mn

    def getArguments(self, m):
        """
        precalculated parameters:
        """
        im = interpolate(m, Function(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % ( str(im)))

        isigma_0 = self.getSigma0(im[0])
        isigma_0_stations =  self.getSigma0(im_stations[:,0])
        iMn = self.getSigma0(im[1])
        
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return im, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, m, im, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)

        gm=grad(m, where=im.getFunctionSpace())
        R = self.w1 / 2 * integrate(length(gm) ** 2)

        if self.theta>0:
            gm0 = gm[0]
            gm1 = gm[1]
            lgm0_2 = length(gm0)**2 + self.m_epsilon**2
            lgm1_2 = length(gm1)**2 + self.m_epsilon**2
            MX = integrate(self.theta * ((lgm0_2 * lgm1_2) - inner(gm0, gm1) ** 2) * (1 / lgm0_2 + 1 / lgm1_2))/2.
        else:
            MX=0
        V = MX + R + misfit_DC + misfit_IP
        self.logger.debug(
                f'misfit ERT, IP; reg ; Xgrad, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e}; {MX:e} = {V:e}')
        self.logger.debug(
                f'ratios ERT, IP; reg ; Xgrad [%] \t=  {misfit_DC/V*100:g}, {misfit_IP/V*100:g};  {R/V*100:g}; {MX/V*100:g}')
        return V

    def getGradient(self, m,  im, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gm=grad(m, where=im.getFunctionSpace())
        X = self.w1 * gm
        DMisfitDsigma_0, DMisfitDMn  = self.getDMisfit(isigma_0, iMn,*args2)

        Dsigma_0Dm0 = self.getDsigma0Dm(isigma_0, im[0])
        DMnDm1 = self.getDMnDm(iMn, im[1])
        Y=Data(0., (2,), im.getFunctionSpace())
        Y[0] = DMisfitDsigma_0 * Dsigma_0Dm0
        Y[1] = DMisfitDMn * DMnDm1

        if self.theta>0:
            gm0 = gm[0]
            gm1 = gm[1]
            lgm0_2 = length(gm0)**2 + self.m_epsilon**2
            lgm1_2 = length(gm1)**2 + self.m_epsilon**2

            X01 = inner(gm0, gm1)
            f = X01 * (1 / lgm0_2   + 1 / lgm1_2 )
            X[0, :] += self.theta * ((1 + X01 ** 2 / lgm0_2 ** 2) * gm0 - f * gm1)
            X[1, :] += self.theta * ((1 + X01 ** 2 / lgm1_2 ** 2) * gm1 - f * gm0)

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, m,  isigma_0,isigma_0_stations,
                                       iMn, *args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        if initializeHessian and self.theta> 0:
            # as restart is disabaled this will not be called:
            gm = grad(m)
            A = self.Hpde.createCoefficient('A')
            gm = args[1]
            gm0 = gm[0]
            gm1 = gm[1]
            lgm0 = length(gm0) + self.m_epsilon
            lgm1 = length(gm1) + self.m_epsilon
            d01 = inner(gm0, gm1)
            f = 1 / lgm0 ** 2 + 1 / lgm1 ** 2
            A[0, :, 0, :] = (self.w1 + self.theta * (1 + d01 ** 2 / lgm0 ** 4)) * kronecker(3) + \
                    self.theta * (2 * d01 / lgm0 ** 4 * (outer(gm0, gm1) + outer(gm1, gm0)) -\
                        4 * d01 ** 2 / lgm0 ** 6 * outer(gm0, gm0) - f * outer(gm1,gm1))
            A[1, :, 1, :] = ((self.w1 + self.theta * (1 + d01 ** 2 / lgm1 ** 4)) * kronecker(3) +\
                        self.theta * (2 * d01 / lgm1 ** 4 * (outer(gm0, gm1) + outer(gm1, gm0)) -\
                                      4 * d01 ** 2 / lgm1 ** 6 * outer(gm1, gm1) - f * outer(gm0,gm0)))
            H01 = self.theta * (2 * d01 * (1 / lgm0 ** 4 * outer(gm0, gm0) + 1 / lgm1 ** 4 * outer(gm1, gm1)) -\
                                f * (d01 * kronecker(3) + outer(gm1, gm0)))
            A[0, :, 1, :] = H01
            A[1, :, 0, :] = transpose(H01)
            self.Hpde.setValue(A=A)

        dm=Data(0., (2,), Solution(self.domain))
        for i in [0,1]:
            self.Hpde.setValue(X=r[1][i], Y=r[0][i], y=r[2][i])
            dm[i] = self.Hpde.getSolution()
            self.logger.debug(f"search direction component {i} = {str(dm[i])}.")
        return dm

    def getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0], m) + inner(r[1], grad(m)))

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)

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