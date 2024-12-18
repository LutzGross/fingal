"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException

from .tools import setupERTPDE, getSourcePotentials, makeMaskForOuterSurface, getAdditivePotentials, DataMisfitQuad, DataMisfitLog, setupPDESystem
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple

class ERTMisfitCostFunction(CostFunction):
    """
    this the data misfit cost function for ERT
    """

    def __init__(self, domain=None, data=None, sigma_src=1., pde_tol=1.e-8,
                 maskZeroPotential = None, stationsFMT="e%s", useLogMisfitDC=False,
                 dataRTolDC=1e-4, logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`. Resistance and error are used.
        :param sigma_src: electric conductivity used to calculate the source potential.
        :param pde_tol: tolerance for the forward and adjoint PDE
        :param maskZeroPotential: mask of where potential is set to zero
        :param useLogMisfitDC: if True, logarithm of data is used in the misfit
        :param dataRTolDC: drop tolerance for small data, that is for data with values smaller than
                            dataRTolDC * max(resistance data).
        :stationsFMT: format string to map station keys A to mesh tags stationsFMT%A or None
        :param logger: `logging.Logger` if None, Logger `fingal` is used/created.
        """
        if logger is None:
            self.logger = logging.getLogger('fingal')
        else:
            self.logger = logger

        super().__init__(**kargs)
        self.domain = domain
        self.stationsFMT = stationsFMT
        if not maskZeroPotential or not Lsup(maskZeroPotential) > 0:
            raise ValueError("non-zero mask for zero potential must be given.")
        self.maskZeroPotential = maskZeroPotential
        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.forward_pde = setupERTPDE(domain)
        self.forward_pde.getSolverOptions().setTolerance(pde_tol)
        self.forward_pde.setValue(q=self.maskZeroPotential)
        self.data = data
        self.useLogMisfitDC=useLogMisfitDC
        self.dataRTolDC = dataRTolDC
        # extract values  by the Locator
        station_locations = [data.getStationLocationByKey(k) for k in data.getStationNumeration() ]
        #self.__grab_values_stations = Locator(Solution(domain), station_locations) # Dirac?
        #TODO: is this really working?
        self.__grab_values_stations = Locator(DiracDeltaFunctions(domain), station_locations)  # Dirac?
        self.sigma_src = sigma_src # needs to be a constant.
        self.setSourcePotentials()  # S_s is in the notes
        self.have_all_adjoints = set(self.data.getObservationStations()).issubset(self.data.getInjectionStations())
        if self.have_all_adjoints:
            self.logger.info("Have all adjoint potentials.")

        # build the misfit data (indexed by source index):
        self.misfit_DC = {}  # potential increment to injection field to get DC potential
        data_atol_DC= self.dataRTolDC * data.getMaximumResistence()

        nd_DC = 0 # counting number of data
        n_small_DC = 0 # number of dropped observations
        for A, B in self.data.injectionIterator():
            iA = self.data.getStationNumber(A)
            iB = self.data.getStationNumber(B)
            obs = self.data.getObservations(A, B)
            # ........... DC part ................................................................
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
        self.logger.info(f"Absolute data cut-off tolerance is {data_atol_DC:e}.")
        self.logger.info(f"{nd_DC} DC data are records used. {n_small_DC} small values found.")
        if nd_DC > 0:
            for iA, iB in self.misfit_DC:
                self.misfit_DC[(iA, iB)].rescaleWeight(1. / nd_DC)
        else:
            raise ValueError("No data for the inversion.")
        self.ignoreERTMisfit(flag=False)

    def ignoreERTMisfit(self, flag=False):
        """
        Switch on/off the misfit in the cost function. This is used for testing.
        """
        self.ignore_ERTmisfit = flag
        if self.ignore_ERTmisfit:
            self.logger.info(f"WARNING: **** Misfit is ignored in cost function ****")

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
        self.source_potential = getSourcePotentials(self.domain, self.sigma_src, self.data, maskZeroPotential =self.maskZeroPotential,
                                                    stationsFMT=self.stationsFMT, logger=self.logger)
        self.source_potentials_stations = { iA : self.grabValuesAtStations(self.source_potential[iA])
                                               for iA in self.source_potential }
    def fitSigmaRef(self):
        """
        finds a new sigma_ref that gives a better data fit then
        """
        if len(self.misfit_DC) == 0:
            raise ValueError("No Data Available")
        beta = 1.
        a = 0.
        b = 0.
        if self.useLogMisfitDC:
            for iA, iB in self.misfit_DC:
                u = self.source_potentials_stations[iA] - self.source_potentials_stations[iB]
                du = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * (self.misfit_DC[(iA, iB)].data - log(abs(du)) ))
                b+= sum(  self.misfit_DC[(iA, iB)].weightings )
            if b>0:
                beta=exp(a/b)
        else:
            for iA, iB in self.misfit_DC:
                u = self.source_potentials_stations[iA] - self.source_potentials_stations[iB]
                du = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * self.misfit_DC[(iA, iB)].data * du )
                b+= sum(  self.misfit_DC[(iA, iB)].weightings * du**2 )
            if b > 0:
                beta = a/b

        assert beta > 0, "conductivity correction factor must be positive."
        return self.sigma_src/beta

    def getERTModelAndResponse(self, sigma_0, sigma_0_stations):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and sigma_0_face
        they should be given on integration nodes.
        """
        self.logger.info("sigma_0 = %s" % (str(sigma_0)))
        self.logger.info("sigma_0_stations = %s - %s " % (min(sigma_0_stations), max(sigma_0_stations)))
        additive_potentials_DC, src_potential_scale= getAdditivePotentials(self.forward_pde,
                                                       sigma = sigma_0,
                                                       schedule = self.data,
                                                       sigma_stations = sigma_0_stations,
                                                       source_potential = self.source_potential,
                                                       sigma_src = self.sigma_src,
                                                       logger=self.logger)
        additive_potentials_DC_stations = {iA : self.grabValuesAtStations(additive_potentials_DC[iA])
                                              for iA in additive_potentials_DC}
        #if self.logger.isEnabledFor(logging.DEBUG) :
        #    for iA in self.source_potential:
        #        self.logger.debug("DC secondary potential %d :%s" % (iA, str(additive_potentials_DC[iA])))
        return additive_potentials_DC, additive_potentials_DC_stations, src_potential_scale

    def getMisfit(self, sigma_0, sigma_0_stations, additive_potentials_DC, additive_potentials_DC_stations, src_potential_scale, *args):
        """
        return the misfit in potential and chargeaiBlity weighted by weightingMisfitDC factor
        """
        potentials_DC_stations = {}
        for iA in additive_potentials_DC_stations:
            potentials_DC_stations[iA] = src_potential_scale[iA] * self.source_potentials_stations[iA] + additive_potentials_DC_stations[iA]
        misfit_0 = 0.
        if not self.ignore_ERTmisfit:
            for iA, iB in self.misfit_DC:
                misfit_0 += self.misfit_DC[(iA, iB)].getValue(
                    potentials_DC_stations[iA] - potentials_DC_stations[iB])
        return misfit_0
    def getDMisfit(self, sigma_0, sigma_0_stations, additive_potentials_DC, additive_potentials_DC_stations, src_potential_scale, *args):
        """
        returns the derivative of the misfit function with respect to sigma_0 with respect to Mn
        """
        DMisfitDsigma_0 = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))
        if not self.ignore_ERTmisfit:
            SOURCES=np.zeros( (self.data.getNumStations(), self.data.getNumStations()), float)
            potentials_DC_stations = {}
            for iA in additive_potentials_DC_stations:
                potentials_DC_stations[iA] = src_potential_scale[iA] * self.source_potentials_stations[iA] + additive_potentials_DC_stations[iA]

            for iA, iB in self.misfit_DC:
                dmis_0 = self.misfit_DC[(iA, iB)].getDerivative(
                    potentials_DC_stations[iA] - potentials_DC_stations[iB])
                for i in range(len(self.misfit_DC[(iA, iB)])):
                    iM = self.misfit_DC[(iA, iB)].iMs[i]
                    iN = self.misfit_DC[(iA, iB)].iNs[i]
                    #ABMN
                    SOURCES[iA, iM] += dmis_0[i]
                    SOURCES[iB, iM] -= dmis_0[i]
                    SOURCES[iA, iN] -= dmis_0[i]
                    SOURCES[iB, iN] += dmis_0[i]
            if self.have_all_adjoints:
                VA = { iA : src_potential_scale[iA] * self.source_potential[iA] + additive_potentials_DC[iA] for iA in additive_potentials_DC.keys()}
                for iA in additive_potentials_DC.keys():
                    RA_star = Scalar(0., Solution(self.domain))
                    for iM in VA:
                        RA_star+= SOURCES[iA, iM] * VA[iM]
                    # self.logger.debug("DC adjoint potential %d :%s" % (iA, str(VA_star)))
                    DMisfitDsigma_0 -= inner(grad(RA_star), grad(VA[iA]))
            else:
                self.forward_pde.setValue(A=sigma_0 * kronecker(self.forward_pde.getDim()), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
                for iA in additive_potentials_DC.keys():
                        s = Scalar(0., DiracDeltaFunctions(self.forward_pde.getDomain()))
                        for M in self.data.getStationNumeration():
                            iM = self.data.getStationNumber(M)
                            if self.stationsFMT is None:
                                s.setTaggedValue(M, SOURCES[iA, iM])
                            else:
                                s.setTaggedValue(self.stationsFMT % M,  SOURCES[iA, iM])
                        self.forward_pde.setValue(y_dirac=s)
                        RA_star=self.forward_pde.getSolution()
                        #self.logger.debug("DC adjoint potential %d :%s" % (iA, str(VA_star)))
                        VA= src_potential_scale[iA] * self.source_potential[iA] + additive_potentials_DC[iA]
                        DMisfitDsigma_0 -= inner(grad(RA_star), grad(VA))
            self.logger.debug("%s adjoint potentials calculated." % len(additive_potentials_DC.keys()))
        return DMisfitDsigma_0


class ERTInversionH1(ERTMisfitCostFunction):
    """
    ERT inversion with H1 regularization
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, m_ref = None, sigma_src=None, w1=1., useL1Norm=False, epsilonL1Norm=1e-4,
                 maskFixedProperty = None, maskZeroPotential = None, fixTop=False,
                 pde_tol=1.e-8, reg_tol=None, stationsFMT="e%s", logclip=5, dataRTolDC=1e-4,
                 useLogMisfitDC=False, logger=None, EPSILON=1e-15, **kargs):
        """
        :param domain: pde domain
        :param data: survey data, is `fingal.SurveyData`
        :param sigma_0_ref: sigma = sigma_0_ref * exp(m)
        :param sigma_src: conductivity used to calculate the injection/source potentials (need to be constant)
        :param w1: weighting H1 regularization  int grad(m)^2
        :param useL1Norm: use L1+eps regularization
        :param epsilonL1Norm: cut-off for L1+eps regularization
        :param maskFixedProperty: mask of region where sigma_0 is fixed to sigma_0_ref * exp(m_ref).
            If not set the all faces are fixed except for the top surface if fixTop is set.
        :param fixTop: sigma_0 is also fixed at the top  if maskFixedProperty is not set.
        :param maskZeroPotential: mask where potentials are set to zero
        :param pde_tol: tolerance for the PDE solver
        :param reg_tol: tolerance for the PDE solver in the regularization
        :param stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :param logclip: cliping for p to avoid overflow in conductivity calculation
        :param EPSILON: absolute tolerance for zero values.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential= maskZeroPotential,
                         stationsFMT=stationsFMT, logger=logger, useLogMisfitDC=useLogMisfitDC,
                         dataRTolDC=dataRTolDC, **kargs)
        self.logclip = logclip
        self.m_ref = m_ref
        # its is assumed that sigma_ref on the faces is not updated!!!
        if maskFixedProperty is None:
            x = self.forward_pde.getDomain().getX()
            qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
            qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
            qz = whereZero(x[2] - inf(x[2]))
            if fixTop:
                qz += whereZero(x[2] - sup(x[2]))
            self.mask_fixed_property  =  qx + qy + qz
        else:
            self.mask_fixed_property = wherePositive(maskFixedProperty)
        self.useL1Norm = useL1Norm
        self.epsilonL1Norm = epsilonL1Norm
        assert not self.useL1Norm, "L1-regularization is not fully implemented yet."
        # regularization

        self.factor2D=1
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.setValue(q=self.mask_fixed_property)
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.setW1(w1)

        #  reference conductivity:
        self.updateSigma0Ref(sigma_0_ref)

    def setW1(self, w1=1):
        self.w1 = w1
        K=kronecker(3)
        K[1,1]*=self.factor2D
        self.Hpde.setValue(A=self.w1 * K)
        self.logger.debug(f'w1 = {self.w1:g}')

    def updateSigma0Ref(self, sigma_0_ref):
        """
        set a new reference conductivity
        """
        self.sigma_0_ref=sigma_0_ref
        self.logger.info("Reference conductivity sigma_0_ref is set to %s" % (str(self.sigma_0_ref)))


    def getSigma0(self, m, applyInterploation=False):
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        sigma_0 = self.sigma_0_ref * exp(im)
        return sigma_0

    def getDsigma0Dm(self, sigma_0, m):
        return sigma_0


    def extractPropertyFunction(self, m):
        if self.m_ref:
            m = dm + self.m_ref
        else:
            m = dm
        return m

    def getArguments(self, dm):
        """
        precalculated parameters:
        """
        m=self.extractPropertyFunction(dm)
        im = interpolate(m, Function(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % ( str(im)))

        isigma_0 = self.getSigma0(im)
        isigma_0_stations =  self.getSigma0(im_stations)
        args2 = self.getERTModelAndResponse(isigma_0, isigma_0_stations)
        return im, isigma_0, isigma_0_stations, args2

    def getValue(self, dm, im, isigma_0, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_stations, *args2)

        gdm=grad(dm, where=im.getFunctionSpace())
        lgdm2=gdm[0]**2+self.factor2D*gdm[1]**2+gdm[2]**2
        if self.useL1Norm:
            R=  self.w1 * integrate(sqrt(lgdm2 + self.epsilonL1Norm ** 2))
        else:
            R = self.w1 / 2 * integrate(lgdm2)

        V = R + misfit_0

        self.logger.debug(
                f'misfit ERT, reg ; total \t=  {misfit_0:e}, {R:e} = {V:e}')
        self.logger.debug(
                f'ratios ERT, reg  [%] \t=  {misfit_0/V*100:g}, {R/V*100:g}')
        return V

    def getGradient(self, dm, im, isigma_0, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gdm=grad(dm, where=im.getFunctionSpace())
        X =gdm
        X[1]*=self.factor2D
        if self.useL1Norm:
            lgdm2 = gdm[0] ** 2 + self.factor2D * gdm[1] ** 2 + gdm[2] ** 2
            X *= 1/sqrt(lgdm2 + self.epsilonL1Norm ** 2)
        X*=self.w1
        DMisfitDsigma_0 = self.getDMisfit(isigma_0, isigma_0_stations, *args2)
        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)

        Y = DMisfitDsigma_0 * Dsigma_0Dm
        return ArithmeticTuple(Y, X )

    def getInverseHessianApproximation(self, r, dm, im, isigma_0, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        clip_property_function=10
        if self.useL1Norm:
            if self.HpdeUpdateCount > 1:
                gdm=grad(dm)
                lgdm2 = gdm[0] ** 2 + self.factor2D * gdm[1] ** 2 + gdm[2] ** 2
                L = sqrt(lgdm2 + self.epsilonL1Norm ** 2)
                self.Hpde.setValue(A=self.w1 * ( 1/L  * kronecker(3) - 1/L**3 * outer(gdm, gdm)) )
                self.HpdeUpdateCount+=1
                print("TO DO")

        self.Hpde.setValue(X=r[1], Y=r[0])
        p = self.Hpde.getSolution()
        self.logger.debug(f"search direction component = {str(p)}.")
        return p
    def getDualProduct(self, dm, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(r[0] * dm + inner(r[1], grad(dm)) )

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
# TODO
class ERTInversionH2(ERTMisfitCostFunction):
    """
    ERT inversion with H2 regularization
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1., m_ref=None,
                 maskZeroPotential = None, fixTop=False, save_memory = False,
                 pde_tol=1.e-8, reg_tol=None, stationsFMT="e%s", logclip=5, dataRTolDC=1e-4,
                 useLogMisfitDC=False, logger=None, EPSILON=1e-15, **kargs):
        """
        :param domain: pde domain
        :param data: survey data, is `fingal.SurveyData`
        :param sigma_0_ref: sigma = sigma_0_ref * exp(m)
        :param sigma_src: conductivity used to calculate the injection/source potentials (need to be constant)
        :param w1: weighting H1 regularization  int grad(m)^2
        :param useL1Norm: use L1+eps regularization
        :param epsilonL1Norm: cut-off for L1+eps regularization
        :param maskZeroPotential: mask of the faces where potential field 'radiation' is set.
        :param pde_tol: tolerance for the PDE solver
        :param reg_tol: tolerance for the PDE solver in the regularization
        :param stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :param logclip: cliping for p to avoid overflow in conductivity calculation
        :param EPSILON: absolute tolerance for zero values.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential= maskZeroPotential, dataRTolDC=dataRTolDC,
                         stationsFMT=stationsFMT, logger=logger, useLogMisfitDC=useLogMisfitDC, **kargs)
        self.logclip = logclip
        self.m_ref = m_ref
        self.sigma_0_ref = sigma_0_ref
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))
        # masks for surfaces:
        x=self.forward_pde.getDomain().getX()
        qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
        qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
        qz = whereZero(x[2] - inf(x[2]))
        if fixTop:
            qz += whereZero(x[2] - sup(x[2]))


        # recovery of propery function from M:
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
        else:
            self.Hpdes = []
            for k, q in enumerate( [ qy + qz, qx + qz,  qx + qy ]):
                pde =  setupERTPDE(self.domain)
                pde.getSolverOptions().setTolerance(reg_tol)
                pde.setValue(q=q)
                self.Hpdes.append(pde)

        self.setW1(w1)
        #  reference conductivity:
        self.updateSigma0Ref(sigma_0_ref)

    def setW1(self, w1):
        if self.save_memory:
            self.Hpde.setValue(A=w1 * kronecker(3))
        else:
            for pde in self.Hpdes:
                pde.setValue(A=w1 * kronecker(3))
        self.w1 = w1

    def updateSigma0Ref(self, sigma_0_ref):
        """
        set a new reference conductivity
        """
        self.sigma_0_ref=sigma_0_ref

    def getSigma0(self, m, applyInterploation=False):
        if applyInterploation:
            im = interpolate(m, Function(self.domain))
        else:
            im = m
        im = clip(im, minval=-self.logclip, maxval=self.logclip)
        sigma_0 = self.sigma_0_ref * exp(im)
        return sigma_0

    def getDsigma0Dm(self, sigma_0, m):
        return sigma_0

    def extractPropertyFunction(self, dM):
        self.mpde.setValue(X=interpolate(dM, Function(dM.getDomain())), Y=Data(), y=Data())
        dm=self.mpde.getSolution()
        if self.m_ref:
            m=dm+self.m_ref
        else:
            m=dm
        return m
    def getArguments(self, dM):
        """
        precalculated parameters:
        """
        m=self.extractPropertyFunction(dM)

        im = interpolate(m, Function(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("dM = %s" % ( str(dM)))
        self.logger.debug("m = %s" % ( str(im)))
        isigma_0 = self.getSigma0(im)
        isigma_0_stations =  self.getSigma0(im_stations)
        args2 = self.getERTModelAndResponse(isigma_0, isigma_0_stations)
        return m, im,  isigma_0, isigma_0_stations, args2

    def getValue(self, dM, m, im, isigma_0, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_stations, *args2)
        gdM=grad(dM, where=im.getFunctionSpace())
        reg = self.w1 / 2 * integrate(length(gdM)**2)
        V = reg + misfit_0

        self.logger.debug(f'misfit ERT, reg - total \t=  {misfit_0:e}, {reg:e} -  {V:e}.')
        self.logger.debug(f'ratios ERT, reg - total  [%] \t=  {misfit_0/V*100:g}, {reg/V*100:g} - 100.')
        return V

    def getGradient(self, dM, m, im, isigma_0, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        Ystar = Scalar(0.,im.getFunctionSpace() )
        gdM=grad(dM, where=im.getFunctionSpace())
        X = gdM * self.w1

        DMisfitDsigma_0 = self.getDMisfit(isigma_0, isigma_0_stations, *args2)
        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Ystar+= DMisfitDsigma_0 * Dsigma_0Dm
        self.mpde.setValue(X=Data(), Y = Ystar)
        Y=grad(self.mpde.getSolution(), where=X.getFunctionSpace())
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, dM, m, im, isigma_0, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        dM = Vector(0., dM.getFunctionSpace())
        if self.save_memory:
            for k in range(3):
                self.Hpde.setValue(q= self.Hpde_qs[k], X=r[1][k,:], Y=r[0][k])
                dM[k] = self.Hpde.getSolution()
                self.logger.debug(f"search direction component {k} = {str(dM[k])}.")
        else:
            for k in range(3):
                self.Hpdes[k].setValue(X=r[1][k,:], Y=r[0][k])
                dM[k] = self.Hpdes[k].getSolution()
                self.logger.debug(f"search direction component {k} = {str(dM[k])}.")
        return dM

    def getDualProduct(self, M, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0] , M) + inner(r[1], grad(M)) )

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)
class ERTInversionGaussBase(ERTMisfitCostFunction):
    """
    ERT inversion with Gaussian regulization with decoupled preconditioning
    """
    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1., length_scale =1.,
                 penalty_factor=1,
                 maskZeroPotential = None, fixTop=False, m_ref = None, save_memory = False,
                 pde_tol=1.e-8, reg_tol=None, stationsFMT="e%s", logclip=5, dataRTolDC=1e-4,
                 useLogMisfitDC=False, EPSILON=1e-15, logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: survey data, is `fingal.SurveyData`
        :param sigma_0_ref: sigma = sigma_0_ref * exp(m)
        :param sigma_src: conductivity used to calculate the injection/source potentials (need to be constant)
        :param w1: weighting H2 regularization int (m+a^2*del(m)) ^2
        :param length_scale: length scale factor a
        :param maskZeroPotential: mask of the faces where potential field 'radiation' is set.
        :param pde_tol: tolerance for the PDE solver
        :param reg_tol: tolerance for the PDE solver in the regularization
        :param stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :param logclip: cliping for p to avoid overflow in conductivity calculation
        :param EPSILON: absolute tolerance for zero values.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         maskZeroPotential= maskZeroPotential, dataRTolDC=dataRTolDC,
                         stationsFMT=stationsFMT, logger=logger, useLogMisfitDC=useLogMisfitDC, **kargs)
        self.logclip = logclip
        self.pde_tol = pde_tol
        self.save_memory = save_memory
        self.m_ref = m_ref
        self.sigma_0_ref = sigma_0_ref
        self.length_scale = length_scale
        self.penalty_factor = penalty_factor
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))
        self.logger.info("length scale a = %s" % (str(self.length_scale)))
        # masks for surfaces:
        x=self.forward_pde.getDomain().getX()
        qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
        qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
        qz = whereZero(x[2] - inf(x[2]))
        if fixTop:
            qz += whereZero(x[2] - sup(x[2]))
        self.q = [ 0*(qx+qy+qz), 0*(qy + qz), 0*(qx + qz),  0*(qx + qy) ]
        # regularization
        if reg_tol is None:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.setUpInitialHessian(w1=w1, reg_tol=reg_tol)
        #  reference conductivity:
        self.updateSigma0Ref(sigma_0_ref)

    def updateSigma0Ref(self, sigma_0_ref):
        """
        set a new reference conductivity
        """
        self.sigma_0_ref=sigma_0_ref

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

    def extractPropertyFunction(self, dM):
        if self.m_ref:
            return dM[0] + self.m_ref
        else:
            return dM[0]

    def getArguments(self, dM):
        """
        precalculated parameters:
        """
        if self.m_ref:
            m = dM[0] + self.m_ref
        else:
            m = dM[0]
        im = interpolate(m, Function(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % ( str(im)))

        isigma_0 = self.getSigma0(im)
        isigma_0_stations =  self.getSigma0(im_stations)
        args2 = self.getERTModelAndResponse(isigma_0, isigma_0_stations)
        return im, isigma_0, isigma_0_stations, args2

    def getValue(self, dM, im, isigma_0, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_stations, *args2)
        a=self.length_scale
        b = self.penalty_factor
        graddm =grad(dM[0], where=im.getFunctionSpace())
        divdM = div(dM[1:], where=im.getFunctionSpace())
        curldM = curl(dM[1:], where = im.getFunctionSpace())

        reg = self.w1 * integrate((dM[0]-a*divdM)**2)
        gM =  self.w1 * b**2 * integrate(length(dM[1:]-a*graddm)**2)
        cM =  self.w1 * integrate(length(a * curldM)**2)

        R = (reg + gM + cM )/2
        V = R + misfit_0
        self.logger.debug(
                f'misfit ERT, reg, div, curl; total \t=  {misfit_0:e}, {reg:e}, {gM:e}, {cM:e}; {V:e}')
        self.logger.debug(
                f'ratios ERT, reg  [%] \t=  {misfit_0/V*100:g}, {R/V*100:g}')
        return V

    def getGradient(self, dM, im, isigma_0, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        a=self.length_scale
        b=self.penalty_factor
        graddm =grad(dM[0], where=im.getFunctionSpace())
        divdM = div(dM[1:], where=im.getFunctionSpace())
        curldM = curl(dM[1:], where = im.getFunctionSpace())

        X = Data(0, (4,3), Function(dM.getDomain()))
        Y = Data(0, (4,),  Function(dM.getDomain()))

        Y[1:]+= b**2 *  (dM[1:]- a * graddm)
        X[0,:] = - a * b**2 *  (dM[1:]- a * graddm)

        d=dM[0]-a*divdM
        Y[0]+= d
        for i in range(0, self.domain.getDim()):
                X[1+i,i]= - a * d

        X[3,1] = a**2 * curldM[0]
        X[2,2] = - X[3,1]
        X[1,2] = a**2 * curldM[1]
        X[3,0] = - X[1,2]
        X[2,0] =  a**2 * curldM[2]
        X[1,1] = - X[2,0]

        X*=self.w1
        Y*=self.w1
        #
        DMisfitDsigma_0  = self.getDMisfit(isigma_0, isigma_0_stations, *args2)
        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Y[0] += DMisfitDsigma_0 * Dsigma_0Dm

        return ArithmeticTuple(Y, X)

    def getDualProduct(self, dM, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0],dM) + inner(r[1], grad(dM)))

    def getNorm(self, dM):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(dM)

    def getSqueezeFactor(self, dM, P):
        if self.m_ref:
            m = dM + self.m_ref
        else:
            m=dM
        a=inf((self.logclip-sign(P[0])*m)/(abs(P[0])+1e-99) )
        if a > 0:
            return a
        else:
            return None


class ERTInversionGaussWithDiagonalHessian(ERTInversionGaussBase):
    def setUpInitialHessian(self, w1 = 1,  reg_tol=None):
        if not reg_tol:
            reg_tol = min(sqrt(self.pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        if self.save_memory:
            self.Hpde = setupERTPDE(self.domain)
            self.Hpde.getSolverOptions().setTolerance(reg_tol)
        else:
            self.Hpdes = []
            for i in range(4):
                pde =  setupERTPDE(self.domain)
                pde.setValue( q= self.q[i] )
                self.Hpdes.append(pde)

        self.setW1(w1)

    def setW1(self, w1):
        a = self.length_scale
        b = self.penalty_factor
        if self.save_memory:
            self.Hpde.setValue(A=w1 * a ** 2 * kronecker(3), D=w1)
        else:
            [ pde.setValue(A=w1 * a ** 2 * kronecker(3), D=w1) for pde in self.Hpdes ]
        self.w1=w1

    def getInverseHessianApproximation(self, r, dM, im, isigma_0, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        w1=self.w1
        a = self.length_scale
        b = self.penalty_factor

        P=Data(0., dM.getShape(), dM.getFunctionSpace() )
        for i in range(dM.getDomain().getDim()+1):
            if self.save_memory:
                if i == 0:
                    self.Hpde.setValue(X=r[1][i,:], Y=r[0][i], q=self.q[i], A=(a**2*b**2*w1)* kronecker(3) , D=w1 )
                else:
                    self.Hpde.setValue(X=r[1][i,:], Y=r[0][i], q=self.q[i], A=(a**2*w1)* kronecker(3) , D=w1*b**2)
                P[i] = self.Hpde.getSolution()
            else:
                self.Hpdes[i].setValue(X=r[1][i,:], Y=r[0][i])
                P[i] = self.Hpdes[i].getSolution()
            self.logger.debug(f"search direction component {i} = {str(P[i])}.")
        return P

class ERTInversionGauss(ERTInversionGaussBase):
    def setUpInitialHessian(self, w1=1, reg_tol=None):
        self.Hpde = setupPDESystem(self.domain, numEquations=4, symmetric=True, tolerance=reg_tol)
        if not reg_tol:
            reg_tol=min(sqrt(self.pde_tol), 1e-3)
        self.Hpde.setValue(q=self.q[0]*[1,0,0,0]+self.q[1]*[0,1,0,0]+self.q[2]*[0,0,1,0]+self.q[3]*[0,0,0,1])
        self.setW1(w1)

    def setW1(self, w1):
        a = self.length_scale
        b = self.penalty_factor
        A = self.Hpde.createCoefficient('A')
        B = self.Hpde.createCoefficient('B')
        C = self.Hpde.createCoefficient('C')
        D = self.Hpde.createCoefficient('D')
        A[0, 0, 0, 0] = a ** 2 * b ** 2
        A[0, 1, 0, 1] = a ** 2 * b ** 2
        A[0, 2, 0, 2] = a ** 2 * b ** 2
        A[1, 0, 1, 0] = a ** 2
        A[1, 0, 2, 1] = a ** 2
        A[1, 0, 3, 2] = a ** 2
        A[1, 1, 1, 1] = a ** 2
        A[1, 1, 2, 0] = -a ** 2
        A[1, 2, 1, 2] = a ** 2
        A[1, 2, 3, 0] = -a ** 2
        A[2, 0, 1, 1] = -a ** 2
        A[2, 0, 2, 0] = a ** 2
        A[2, 1, 1, 0] = a ** 2
        A[2, 1, 2, 1] = a ** 2
        A[2, 1, 3, 2] = a ** 2
        A[2, 2, 2, 2] = a ** 2
        A[2, 2, 3, 1] = -a ** 2
        A[3, 0, 1, 2] = -a ** 2
        A[3, 0, 3, 0] = a ** 2
        A[3, 1, 2, 2] = -a ** 2
        A[3, 1, 3, 1] = a ** 2
        A[3, 2, 1, 0] = a ** 2
        A[3, 2, 2, 1] = a ** 2
        A[3, 2, 3, 2] = a ** 2
        B[0, 0, 1] = -a * b ** 2
        B[0, 1, 2] = -a * b ** 2
        B[0, 2, 3] = -a * b ** 2
        B[1, 0, 0] = -a
        B[2, 1, 0] = -a
        B[3, 2, 0] = -a
        C[0, 1, 0] = -a
        C[0, 2, 1] = -a
        C[0, 3, 2] = -a
        C[1, 0, 0] = -a * b ** 2
        C[2, 0, 1] = -a * b ** 2
        C[3, 0, 2] = -a * b ** 2
        D[0, 0] = 1
        D[1, 1] = b ** 2
        D[2, 2] = b ** 2
        D[3, 3] = b ** 2
        self.Hpde.setValue(A= w1 * A, B = w1 *B, C = w1 * C, D = w1 *D)
        self.w1=w1

    def getInverseHessianApproximation(self, r, dM, im, isigma_0, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(X=r[1], Y=r[0])
        P = self.Hpde.getSolution()
        for i in range(dM.getDomain().getDim()+1):
            self.logger.debug(f"search direction component {i} = {str(P[i])}.")
        return P


