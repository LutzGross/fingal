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
            self.logger = logging.getLogger('fingal.IPInversion')
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
        error_DC_max, error_IP_max = 0., 0.
        error_DC2_invsum, error_IP2_invsum = 0., 0.
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

                error_DC_max = max( error_DC_max, max(error_DC))
                error_DC2_invsum+= sum(abs(data_DC)**2/error_DC**2)
                print("DC:", data_DC, error_DC, data_DC/error_DC, data.getMaximumResistence()/error_DC, error_DC_max, error_DC2_invsum)
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
                    #C = np.array([self.data.getChargeabilityData((A, B, M, N)) for M, N in obs_IP])
                    error_IP = np.array([self.data.getSecondaryResistenceError((A, B, M, N)) for M, N in obs_IP]) + data_atol_IP
                    self.misfit_IP[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_IP, iNs=iNs, injections=(A, B),
                                                              weightings=1. / error_IP ** 2 / n_use_IP)
                error_IP_max = max( error_IP_max, max(error_IP))
                error_IP2_invsum+= sum(data_IP**2/error_IP**2)
                #print("IP:", data_IP, error_IP,  data_IP/error_IP, data.getMaximumSecondaryResistence()/error_IP, error_IP_max, error_IP2_invsum)
                #print("R:", data_IP/data_DC, C, data_IP/error_IP *  error_DC/data_DC )
                nd_IP += len(self.misfit_IP[(iA, iB)])
        print("f_DC = ",error_DC2_invsum/nd_DC)
        print("f_IP = ", error_IP2_invsum/nd_IP)
        print(1/data.getMaximumChargeability()**2, error_IP2_invsum/error_DC2_invsum
              )
        self.logger.info(f"Data drop tolerance for resistance is {data_atol_DC:e}.")
        self.logger.info(f"{nd_DC} DC data are records used. {n_small_DC} small values found.")
        self.logger.info(f"Data drop tolerance for secondary resistance is {data_atol_IP:e}.")
        self.logger.info(f"{nd_IP} IP data are records used. {n_small_IP} small values found.")
        #ffff2 = data.getMaximumChargeability() ** 2
        ffff2 = 1
        if not nd_IP + nd_DC >0 :
            raise ValueError("No data for the inversion.")
        if nd_DC > 0:
            for iA, iB in self.misfit_DC:
                self.misfit_DC[(iA, iB)].rescaleWeight(1. / (2 * nd_DC))
        if nd_IP > 0:
            for iA, iB in self.misfit_IP:
                self.misfit_IP[(iA, iB)].rescaleWeight(1. / (2 * nd_IP)  * ffff2  )
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

    def fitSigmaAndMnRef(self):
        """
        finds a new sigma_ref that gives a better data fit then
        """
        if len(self.misfit_DC) == 0:
            raise ValueError("No Data Available")
        betas = []
        for dataset, uselog in [(self.misfit_DC, self.useLogMisfitDC), (self.misfit_IP, self.useLogMisfitIP)]:
            beta = 1.
            a, b= 0.,0.
            if uselog:
                for iA, iB in dataset:
                    u = self.source_potentials_stations[iA] - self.source_potentials_stations[iB]
                    du = dataset[(iA, iB)].getDifference(u)
                    a+= sum(  dataset[(iA, iB)].weightings * (dataset[(iA, iB)].data - log(abs(du)) ))
                    b+= sum(  dataset[(iA, iB)].weightings )
                if b>0:
                    beta=exp(a/b)
            else:
                for iA, iB in dataset:
                    u = self.source_potentials_stations[iA] - self.source_potentials_stations[iB]
                    du = dataset[(iA, iB)].getDifference(u)
                    a+= sum(  dataset[(iA, iB)].weightings * dataset[(iA, iB)].data * du )
                    b+= sum(  dataset[(iA, iB)].weightings * du**2 )
                if b > 0:
                    beta = a/b
            betas.append(beta)
        print(betas)
        assert betas[0] > 0, "conductivity correction factor must be positive."
        assert abs(betas[0]-betas[1]) > 0, "divison by zero in chargeability."

        sigma_ref = self.sigma_src/betas[0]
        Mn_ref = betas[1] / (betas[0] - betas[1]) * sigma_ref
        assert sigma_ref >0, f"new sigma_ref={sigma_ref} must be positive."
        assert Mn_ref > 0, f"new Mn_ref={Mn_ref} must be positive."
        return sigma_ref, Mn_ref

    def getIPModelAndResponse(self,sigma_0, sigma_0_stations, Mn):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and chargeaiBlity Mn
        they should be given on integration nodes.
        """
        sigma_oo = sigma_0 + Mn

        self.logger.debug("sigma_0 = %s" % (str(sigma_0)))
        self.logger.debug("sigma_0_stations = %g - %g " % (min(sigma_0_stations), max(sigma_0_stations)))
        self.logger.debug("Mn = %s" % (str(Mn)))
        self.logger.debug("sigma_oo = %s" % (str(sigma_oo)))


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
        if misfit_DC > 0:
            self.logger.info(f"Misfits: DC = {misfit_DC}, IP = {misfit_IP} (ratio = {misfit_IP/misfit_DC})." )
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
            self.DC_pde.setValue(y_dirac=s, X= Mn * grad(WA_star))
            VA_star = self.DC_pde.getSolution()
            UA = VA - WA
            DMisfitDsigma_0 -= inner(grad(VA_star), grad(VA)) + inner(grad(WA_star), grad(WA))
            DMisfitDMn += inner(grad(WA_star), grad(UA))
        # -------------------------
        self.logger.debug("%s adjoint potentials calculated." % len(additive_potentials_DC.keys()))
        return DMisfitDsigma_0,  DMisfitDMn


class IPInversionH1(IPMisfitCostFunction):
    """
   Class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
   using H1-regularization on m-m_ref with m=[m[0], m[1]]:

            m[0]=log(sigma/sigma_ref), m[1]=log(Mn/Mn_ref)

   cross-gradient is added to cost function with weighting factor theta (=0 by default)
   """

    def __init__(self, domain, data, maskZeroPotential, sigma_src=None, pde_tol= 1e-8,
                    stationsFMT="e%s", m_ref = None, fix_top = False, zero_mean_m = False,
                    useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                    weightingMisfitDC=1, sigma_0_ref=1e-4, Mn_ref = 1.e-6,
                    w1=1, theta=0., maskFixedProperty = None,
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
        :param useLogMisfitDC: if set logarithm of DC data is used in misfit.
        :param useLogMisfitIP: if set logarithm of secondary potential (IP) data is used in misfit.
        :param dataRTolDC: relative tolerance for damping small DC data on misfit
        :param dataRTolIP: relative tolerance for damping small secondary potential data in misfit
        :param weightingMisfitDC: weighting factor for DC data in misfit. weighting factor for IP data is one.
        :param sigma_0_ref: reference DC conductivity
        :param Mn_ref: reference Mn conductivity
        :param w1: regularization weighting factor(s) (scalar or numpy.ndarray)
        :param theta: weigthing factor x-grad term
        :param maskFixedProperty: mask of regions where m[0]-m_ref[0] and m[1]-m_ref[1] are fixed.
                                    if not set left, right, front, back and bottom are fixed.
        :param logclip: value of m are clipped to stay between -logclip and +logclip
        :param m_epsilon: threshold for small m, grad(m) values.
        :param reg_tol: tolerance for PDE solve for regularization.
        :param logger: the logger, if not set, 'fingal.IPInversion.H1' is used.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        if logger is None:
            self.logger = logging.getLogger('fingal.IPInversion.H1')
        else:
            self.logger = logger
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol= pde_tol,
                            maskZeroPotential=maskZeroPotential, stationsFMT=stationsFMT,
                            useLogMisfitDC=useLogMisfitDC, dataRTolDC=dataRTolDC,
                            useLogMisfitIP=useLogMisfitIP, dataRTolIP=dataRTolIP,
                            weightingMisfitDC=weightingMisfitDC,
                            logger=self.logger, **kargs)

        self.logclip = logclip
        self.m_epsilon = m_epsilon
        self.m_ref = m_ref
        self.zero_mean_m = zero_mean_m
        # its is assumed that sigma_ref and Mn_ref on the faces is not updated!!!
        if maskFixedProperty is None:
            x = self.DC_pde.getDomain().getX()
            qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
            qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
            qz = whereZero(x[2] - inf(x[2]))
            if fix_top:
                qz += whereZero(x[2] - sup(x[2]))
            self.mask_fixed_property  =   qx + qy + qz
        else:
            self.mask_fixed_property = wherePositive(maskFixedProperty)
        self.updateSigma0Ref(sigma_0_ref)
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
        self.setW1andTheta(w1, theta)

    def setW1andTheta(self, w1, theta):
        self.w1=w1 * np.ones((2,))
        self.theta = theta
        self.logger.debug(f'w1[0] = {self.w1[0]:g}')
        self.logger.debug(f'w1[1] = {self.w1[1]:g}')
        # ....
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
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m[0] = %s" % (str(im[0])))
        self.logger.debug("m[1] = %s" % (str(im[1])))

        isigma_0 = self.getSigma0(im[0])
        isigma_0_stations =  self.getSigma0(im_stations[:,0])
        iMn = self.getMn(im[1])
        
        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return im, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, dm, im, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)

        gdm=grad(dm, where=im.getFunctionSpace())
        print("grad m[0] :", str(length(gdm[0])))
        print("grad m[1] :", str(length(gdm[1])))
        R = 1. / 2 * integrate(self.w1[0] * length(gdm[0]) ** 2 + self.w1[1] * length(gdm[1]) ** 2)

        if self.theta>0:
            gdm0 = gdm[0]
            gdm1 = gdm[1]
            lgdm0_2 = length(gdm0)**2 + self.m_epsilon**2
            lgdm1_2 = length(gdm1)**2 + self.m_epsilon**2
            MX = integrate(self.theta * ((lgdm0_2 * lgdm1_2) - inner(gdm0, gdm1) ** 2) * (1 / lgdm0_2 + 1 / lgdm1_2))/2.
        else:
            MX=0
        V = MX + R + misfit_DC + misfit_IP
        self.logger.debug(
                f'misfit ERT, IP; reg ; Xgrad, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e}; {MX:e} = {V:e}')
        self.logger.debug(
                f'ratios ERT, IP; reg ; Xgrad [%] \t=  {misfit_DC/V*100:g}, {misfit_IP/V*100:g};  {R/V*100:g}; {MX/V*100:g}')
        return V

    def getGradient(self, dm,  im, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gdm=grad(dm, where=im.getFunctionSpace())
        X = Data(0., (2,3), gdm.getFunctionSpace() )
        X[0] = self.w1[0] * gdm[0]
        X[1] = self.w1[1] * gdm[1]
        DMisfitDsigma_0, DMisfitDMn  = self.getDMisfit(isigma_0, iMn,*args2)
        print("DMisfitDsigma_0 ", str(DMisfitDsigma_0))
        print("DMisfitDMn ", str(DMisfitDMn))

        Dsigma_0Dm0 = self.getDsigma0Dm(isigma_0, im[0])
        DMnDm1 = self.getDMnDm(iMn, im[1])
        Y=Data(0., (2,), im.getFunctionSpace())
        Y[0] = DMisfitDsigma_0 * Dsigma_0Dm0
        Y[1] = DMisfitDMn * DMnDm1

        if self.theta>0:
            gdm0 = gdm[0]
            gdm1 = gdm[1]
            lgdm0_2 = length(gdm0)**2 + self.m_epsilon**2
            lgdm1_2 = length(gdm1)**2 + self.m_epsilon**2

            X01 = inner(gdm0, gdm1)
            f = X01 * (1 / lgdm0_2   + 1 / lgdm1_2 )
            X[0, :] += self.theta * ((1 + X01 ** 2 / lgdm0_2 ** 2) * gdm0 - f * gdm1)
            X[1, :] += self.theta * ((1 + X01 ** 2 / lgdm1_2 ** 2) * gdm1 - f * gdm0)

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, dm,  isigma_0,isigma_0_stations,
                                       iMn, *args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        if False and initializeHessian and self.theta> 0:
            # as restart is disabaled this will not be called:
            gdm = grad(dm)
            A = self.Hpde.createCoefficient('A')
            gdm = args[1]
            gdm0 = gdm[0]
            gdm1 = gdm[1]
            lgdm0 = length(gdm0) + self.m_epsilon
            lgdm1 = length(gm1) + self.m_epsilon
            d01 = inner(gdm0, gdm1)
            f = 1 / lgdm0 ** 2 + 1 / lgdm1 ** 2
            A[0, :, 0, :] = (self.w1 + self.theta * (1 + d01 ** 2 / lgdm0 ** 4)) * kronecker(3) + \
                    self.theta * (2 * d01 / lgdm0 ** 4 * (outer(gdm0, gdm1) + outer(gdm1, gdm0)) -\
                        4 * d01 ** 2 / lgdm0 ** 6 * outer(gdm0, gdm0) - f * outer(gdm1,gdm1))
            A[1, :, 1, :] = ((self.w1 + self.theta * (1 + d01 ** 2 / lgm1 ** 4)) * kronecker(3) +\
                        self.theta * (2 * d01 / lgdm1 ** 4 * (outer(gdm0, gdm1) + outer(gdm1, gdm0)) -\
                                      4 * d01 ** 2 / lgdm1 ** 6 * outer(gdm1, gdm1) - f * outer(gdm0,gdm0)))
            H01 = self.theta * (2 * d01 * (1 / lgdm0 ** 4 * outer(gdm0, gdm0) + 1 / lgdm1 ** 4 * outer(gdm1, gdm1)) -\
                                f * (d01 * kronecker(3) + outer(gdm1, gdm0)))
            A[0, :, 1, :] = H01
            A[1, :, 0, :] = transpose(H01)
            self.Hpde.setValue(A=A)
        dm=Data(0., (2,), Solution(self.domain))
        if self.zero_mean_m:
            for i in [0,1]:
                dm[i] = self.Hpde.getSolution(X=r[1][i], Y=r[0][i] )/(self.w1[i] + self.theta)
                self.logger.debug(f"search direction component {i} = {str(dm[i])}.")
        else:
            for i in [0,1]:
                self.Hpde.setValue(X=r[1][i], Y=r[0][i])
                dm[i] = self.Hpde.getSolution()/(self.w1[i] + self.theta)
                self.logger.debug(f"search direction component {i} = {str(dm[i])}.")
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


class IPInversionH2(IPMisfitCostFunction):
    """
   Class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
   using H2-regularization on M = grad(m-m_ref) (as 6 components) with m=[m[0], m[1]]:

            m[0]=log(sigma/sigma_ref), m[1]=log(Mn/Mn_ref)

   To recover m from we solve

            (grad(v), grad(m[i]) = (grad(v), M[i*3:(i+1)*3])

   cross-gradient over V is added to cost function with weighting factor theta (=0 by default)
   """

    def __init__(self, domain, data, maskZeroPotential, sigma_src=None, pde_tol=1e-8,
                 stationsFMT="e%s", m_ref=None, length_scale=None,
                 useLogMisfitDC=False, dataRTolDC=1e-4, useLogMisfitIP=False, dataRTolIP=1e-4,
                 weightingMisfitDC=1, sigma_0_ref=1e-4, Mn_ref=1.e-6,
                 w1=1, theta=0., fixTop=False,  zero_mean_m = True,
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


        self.setW1andTheta(w1, theta)

    def setW1andTheta(self, w1, theta):
        self.w1 = w1 * np.ones((2,))
        self.theta = theta
        self.logger.debug(f'w1[0] = {self.w1[0]:g}')
        self.logger.debug(f'w1[1] = {self.w1[1]:g}')
        self.logger.debug(f'theta = {self.theta:g}')

    def updateSigma0Ref(self, sigma_0_ref):
        """
        set a new reference conductivity
        """
        self.sigma_0_ref = sigma_0_ref
        self.logger.info("Reference conductivity sigma_0_ref is set to %s" % (str(self.sigma_0_ref)))

    def updateMnRef(self, Mn_ref):
        """
        set a new reference conductivity
        """
        self.Mn_ref = Mn_ref
        self.logger.info("Reference normalised chargeability Mn_ref = %s" % (str(self.Mn_ref)))

    def getSigma0(self, m, applyInterploation=False):
        if hasattr(m, "getShape") and not m.getShape() == ():
            m = m[0]
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
            m = m[1]
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
        dm = Data(0.,(2,), Solution(self.domain))

        if self.zero_mean_m:
            dm[0] = self.mpde.getSolution(X=interpolate(dM[0:3], Function(dM.getDomain())))
            dm[1] = self.mpde.getSolution(X=interpolate(dM[3:], Function(dM.getDomain())))
        else:
            self.mpde.setValue(X=interpolate(dM[0:3], Function(dM.getDomain())), Y=Data(), y=Data())
            dm[0]=self.mpde.getSolution()
            self.mpde.setValue(X=interpolate(dM[3:6], Function(dM.getDomain())), Y=Data(), y=Data())
            dm[1] = self.mpde.getSolution()
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
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m[0] = %s" % (str(im[0])))
        self.logger.debug("m[1] = %s" % (str(im[1])))

        isigma_0 = self.getSigma0(im[0])
        isigma_0_stations = self.getSigma0(im_stations[:, 0])
        iMn = self.getMn(im[1])

        args2 = self.getIPModelAndResponse(isigma_0, isigma_0_stations, iMn)
        return im, isigma_0, isigma_0_stations, iMn, args2

    def getValue(self, dM, im, isigma_0, isigma_0_stations, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_DC, misfit_IP = self.getMisfit(*args2)

        gdM = grad(dM, where=im.getFunctionSpace())
        idM = interpolate(dM, im.getFunctionSpace() )
        R = 1. / 2. * integrate(self.w1[0] * ( length(gdM[0:3]) ** 2 + self.a * length(idM[0:3])**2 )
                             + self.w1[1] * ( length(gdM[3:] ) ** 2 + self.a * length(idM[3:])**2 ) )

        if self.theta > 0:
            gdm0 = interpolate(dM[0:3], Function(self.domain))
            gdm1 = interpolate(dM[3:], Function(self.domain))
            lgdm0_2 = length(gdm0) ** 2 + self.m_epsilon ** 2
            lgdm1_2 = length(gdm1) ** 2 + self.m_epsilon ** 2
            MX = integrate(
                self.theta * ((lgdm0_2 * lgdm1_2) - inner(gdm0, gdm1) ** 2) * (1 / lgdm0_2 + 1 / lgdm1_2)) / 2.
        else:
            MX = 0
        V = MX + R + misfit_DC + misfit_IP
        if self.theta > 0:
            self.logger.debug(
                f'misfit ERT, IP; reg ; Xgrad, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e}; {MX:e} = {V:e}')
            self.logger.debug(
                f'ratios ERT, IP; reg ; Xgrad [%] \t=  {misfit_DC / V * 100:g}, {misfit_IP / V * 100:g};  {R / V * 100:g}; {MX / V * 100:g}')
        else:
            self.logger.debug(
                f'misfit ERT, IP; reg, total \t=  {misfit_DC:e}, {misfit_IP:e};  {R:e} = {V:e}')
            self.logger.debug(
                f'ratios ERT, IP; reg \t=  {misfit_DC / V * 100:g}, {misfit_IP / V * 100:g};  {R / V * 100:g}')

        return V

    def getGradient(self, dM, im, isigma_0, isigma_0_stations, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gdM = grad(dM, where=im.getFunctionSpace())
        idM = interpolate(dM, im.getFunctionSpace() )
        Y = Data(0., (6,), im.getFunctionSpace())
        X = Data(0., (6, 3), gdM.getFunctionSpace())

        X[0:3] = self.w1[0] * gdM[0:3]
        X[3:]  = self.w1[1] * gdM[3:]
        Y[0:3] = self.w1[0] * self.a * idM[0:3]
        Y[3:]  = self.w1[1] * self.a * idM[3:]

        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)
        Dsigma_0Dm0 = self.getDsigma0Dm(isigma_0, im[0])
        DMnDm1 = self.getDMnDm(iMn, im[1])
        Ystar0= DMisfitDsigma_0 * Dsigma_0Dm0
        Ystar1= DMisfitDMn * DMnDm1
        if self.zero_mean_m:
            dmstar0 = self.mpde.getSolution( Y = Ystar0 )
            dmstar1 = self.mpde.getSolution( Y = Ystar1 )
        else:
            self.mpde.setValue(X=Data(), Y = Ystar0)
            dmstar0=self.mpde.getSolution()
            self.mpde.setValue(X=Data(), Y = Ystar1)
            dmstar1=self.mpde.getSolution()

        Y[0:3]+= grad(dmstar0, where=X.getFunctionSpace())
        Y[3:]+= grad(dmstar1, where=X.getFunctionSpace())

        if self.theta > 0:
            gdm0 = interpolate(dM[0:3], Function(self.domain))
            gdm1 = interpolate(dM[3:], Function(self.domain))
            lgdm0_2 = length(gdm0) ** 2 + self.m_epsilon ** 2
            lgdm1_2 = length(gdm1) ** 2 + self.m_epsilon ** 2

            X01 = inner(gdm0, gdm1)
            f = X01 * (1 / lgdm0_2 + 1 / lgdm1_2)
            Y[0:3] += self.theta * ((1 + X01 ** 2 / lgdm0_2 ** 2) * gdm0 - f * gdm1)
            Y[3:] += self.theta * ((1 + X01 ** 2 / lgdm1_2 ** 2) * gdm1 - f * gdm0)

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, dM, isigma_0, isigma_0_stations,
                                       iMn, *args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        if False and initializeHessian and self.theta > 0: # THIS DOES NOT WORK
            # as restart is disabaled this will not be called:
            gdm = grad(dm)
            A = self.Hpde.createCoefficient('A')
            gdm = args[1]
            gdm0 = gdm[0]
            gdm1 = gdm[1]
            lgdm0 = length(gdm0) + self.m_epsilon
            lgdm1 = length(gm1) + self.m_epsilon
            d01 = inner(gdm0, gdm1)
            f = 1 / lgdm0 ** 2 + 1 / lgdm1 ** 2
            A[0, :, 0, :] = (self.w1 + self.theta * (1 + d01 ** 2 / lgdm0 ** 4)) * kronecker(3) + \
                            self.theta * (2 * d01 / lgdm0 ** 4 * (outer(gdm0, gdm1) + outer(gdm1, gdm0)) - \
                                          4 * d01 ** 2 / lgdm0 ** 6 * outer(gdm0, gdm0) - f * outer(gdm1, gdm1))
            A[1, :, 1, :] = ((self.w1 + self.theta * (1 + d01 ** 2 / lgm1 ** 4)) * kronecker(3) + \
                             self.theta * (2 * d01 / lgdm1 ** 4 * (outer(gdm0, gdm1) + outer(gdm1, gdm0)) - \
                                           4 * d01 ** 2 / lgdm1 ** 6 * outer(gdm1, gdm1) - f * outer(gdm0, gdm0)))
            H01 = self.theta * (2 * d01 * (1 / lgdm0 ** 4 * outer(gdm0, gdm0) + 1 / lgdm1 ** 4 * outer(gdm1, gdm1)) - \
                                f * (d01 * kronecker(3) + outer(gdm1, gdm0)))
            A[0, :, 1, :] = H01
            A[1, :, 0, :] = transpose(H01)
            self.Hpde.setValue(A=A)

        P = Data(0., (6,), Solution(self.domain))
        for k in [0, 1, 2]:
            if self.save_memory:
                self.Hpde.setValue(q=self.Hpde_qs[k])
            for i in [0, 1]:
                if self.save_memory:
                    self.Hpde.setValue(X=r[1][i*3+k], Y=r[0][i*3+k])
                    P[i * 3 + k] = self.Hpde.getSolution() / (self.w1[i] + self.theta)
                else:
                    self.Hpdes[k].setValue(X=r[1][i * 3 + k], Y=r[0][i * 3 + k])
                    P[i*3+k] = self.Hpdes[k].getSolution() / (self.w1[i] + self.theta)
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
