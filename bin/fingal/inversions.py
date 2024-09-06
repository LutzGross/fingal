"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException
from .tools import setupERTPDE, getSourcePotentials, makeMaskForOuterSurface, getSecondaryPotentials
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple

class DataMisfit(object):
    """
     defines a misfit function F as sum of components  f(diff, data, weighting)>=0 with

     f = sum_i f(diff[i], data[i]], weighting[i])

    which is minimised during an inversion. The diff is given has differences of a potential u
    on electrodes M and N.
    """
    EPS=1E-30
    def __init__(self, iMs, data, iNs=None, weightings=None, injections=(), **kwargs):
        """
        :param iMs: index of electrodes
        :param data:data vectopr of observations
        :param iNs: index of second electrodes. if None monopole data are used
        :param weightings: weighting factors; typically 1/(relative error)**2 or zero
        :param injections: injection pair (not used))
        """
        self.iMs = iMs
        self.iNs = iNs
        self.injections = injections
        # assert not ( weightings is None or quadWeight is None )

        if isinstance(data, np.ndarray):
            self.data = data
        else:
            self.data = np.array(data)

        if isinstance(weightings, np.ndarray):
            self.weightings = weightings
        elif weightings is None:
            self.weightings = np.ones(self.data.shape)
        elif isinstance(weightings, float) or isinstance(weightings, int):
            self.weightings = np.zeros(self.data.shape) + weightings
        else:
            self.weightings = np.array(weightings)
    def __len__(self):
        """
        returns number of observations
        """
        return len(self.iMs)

    def rescaleWeight(self, factor):
        self.weightings *= factor
    def getDifference(self, u):
        """
        differences at electrodes from voltage u
        """
        if self.iNs:
            diff = u[self.iMs] - u[self.iNs]
        else:
            diff = u[self.iMs]
        return diff
    def getValue(self, u):
        """
        return the misfit value for potential u
        """
        diff = self.getDifference(u)
        raise NotImplemented
    def getDerivative(self, u):
        """
        return the derivative with respect to vector  diff[i]
        """
        diff = self.getDifference(u)
        raise NotImplemented

class DataMisfitQuad(DataMisfit):
    """
    A mechanism to quadratic handel misfits

    1/2 * sum weighting[i]*(data[i]-diff)**2
    """
    def getValue(self, u):
        res = self.data-self.getDifference(u)
        dd = abs(res) ** 2 * self.weightings
        return 0.5 * sum(dd)

    def getDerivative(self, u):
        res = self.data-self.getDifference(u)
        return -res * self.weightings

class DataMisfitLog(DataMisfit):
    """
    A mechanism to quadratic handel misfits

   1/2 sum weighting[i]*(log(abs(observations[i])-log(abs(data[i])))**2

    where observations are collected from a value vector iMs or difference using iMs and iNs
    """
    def __init__(self, iMs, data, iNs=None, weightings=None, injections=(), **kwargs):
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        data_log = log(abs(data)+self.EPS)
        super().__init__(iMs, data_log, iNs, weightings, injections, **kwargs)
    def getValue(self, u):
        res = self.data-log(abs(self.getDifference(u))+self.EPS)
        dd = abs(res) ** 2 * self.weightings
        return 0.5 * sum(dd)

    def getDerivative(self, u):
        nn=self.getDifference(u)
        res=self.data-log(abs(nn)+self.EPS)
        return res * self.weightings/nn


class IPMisfitCostFunction(CostFunction):
    """
    this the data misfit costfunction for IP
    """

    def __init__(self, domain=None, data=None, sigma_src=1., pde_tol=1.e-8,
                 mask_outer_faces = None, stationsFMT="e%s", useLogMisfitDC=False,
                 useLogMisfitIP=False, weighting_misfit_DC = 0.5,
                 logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`
        :param sigma_src: conductivity to calculate potential of sources
        :param pde_tol:tolerance for the PDE solve
        :param mask_outer_field: mask of surface to apply 'radiation' condition.
        :param weighting_misfit_DC: weighting factor for DC part of cost function
        :param stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :param useLogMisfitDC: use log of data in DC misfit term
        :param useLogMisfitIP: use log of data in IP misfit term
        :logclip: cliping for p to avoid overflow in conductivity calculation
        """
        if logger is None:
            self.logger = logging.getLogger('fingal')
        else:
            self.logger = logger
        self.domain = domain
        self.datatol = 1e-30
        self.stationsFMT = stationsFMT
        self.mask_outer_faces = makeMaskForOuterSurface(domain, facemask=mask_outer_faces)
        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.forward_pde = setupERTPDE(domain)
        self.forward_pde.getSolverOptions().setTolerance(pde_tol)
        # self.forward_pde.getSolverOptions().setTrilinosParameter("reVs_ooe: type","full")
        self.data = data
        self.useLogMisfitERT = useLogMisfitERT
        self.useLogMisfitIP = useLogMisfitIP
        self.setMisfitWeighting(weighting_misfit_DC)

        # extract values  by the Locator
        station_locations = []
        for k in data.getStationNumeration():
            station_locations.append(data.getStationLocationByKey(k))
        self.__grab_values_at_stations = Locator(Solution(domain), station_locations)  # Dirac?
        # self.grabValues=Locator(DiracDeltaFunctions(domain), station_locations)

        self.sigma_src = sigma_src  # needs to be a constant.
        self.setSourcePotentials()  # S_s is in the notes
        # build the misfit data (indexed by source index):
        self.misfit_IP = {}  # secondary potential
        self.misfit_DC = {}  # potential increment to injection field to get DC potential
        nd0 = 0  # counting number of data
        nd2 = 0
        for A, B in self.data.injectionIterator():
            data_DC = np.array([self.data.getResistenceData((A, B, M, N)) for M, N in obs])
            if self.useLogMisfitERT:
                rel_error_DC = np.array([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs])
                self.misfit_DC[(iA, iB)] = DataMisfitLog(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B),
                                             weightings=1. / rel_error_DC ** 2)
            else:
                self.misfit_DC[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B),
                                            weightings = abs(data_DC) ** 2)

            nd0 += len(self.misfit_DC[(iA, iB)])

            data_IP = []
            w_IP = []
            for M, N in obs:
                v = self.data.getSecondaryResistenceData((A, B, M, N))  # u_0-u_oo
                if self.data.isUndefined(v):
                    data_IP.append(-1000.)
                    w_IP.append(0.)
                else:
                    data_IP.append(v)
                    if self.useLogMisfitERT:
                        e = self.data.getSecondaryResistenceRelError((A, B, M, N))
                        w_IP.append((1 / e) ** 2)
                    else:
                        w_IP.append(abs(v)**2)
                    nd2+=1
            data_IP, w_IP = np.array(data_IP), np.array(w_IP)
            if self.useLogMisfitERT:
                self.misfit_DC[(iA, iB)] = DataMisfitLog(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B), weightings = w_IP)
            else:
                self.misfit_DC[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B), weightings = w_IP)

        self.logger.info("%d/%d DC/IP data records detected." % (nd0, nd2))
        nd = nd0 + nd2
        if nd > 0:
            for iA, iB in self.misfit_IP:
                self.misfit_DC[(iA, iB)].rescaleWeight(1. / nd)
                self.misfit_IP[(iA, iB)].rescaleWeight(1. / nd)

    # DOTO:
    def setMisfitWeighting(self, weighting_misfit_ERT=0):

        self.weightingMisfit_2 = max(0, 1 - weighting_misfit_ERT)
        self.weightingMisfit_0 = max(0, weighting_misfit_ERT)

        #self.weightingMisfit_2=0
        #self.weightingMisfit_0=1

        # assert self.weightingMisfit_2>0 or self.weightingMisfit_0>0
    def grabValuesAtStations(self, m):
        """
        returns the values of m at the stations. The order is given by data.getStationNumeration()
        """
        return np.array(self.__grab_values_at_stations(m))
    def setSourcePotentials(self):
        """
        return the primary electric potential for all injections (A,B) Vs_ooing conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_Field
        """
        self.logger.debug("source conductivity sigma_src =  %s" % (str(self.sigma_src)))
        self.source_potential = getSourcePotentials(self.domain, self.sigma_src, self.data, mask_outer_faces =self.mask_outer_faces, stationsFMT=self.stationsFMT)
        self.source_potentials_at_stations = { iA : self.grabValuesAtStations(self.source_potential[iA])
                                               for iA in self.source_potential }

    def getIPModelAndResponse(self,sigma_0, sigma_0_faces, sigma_0_at_stations, Mn, Mn_faces, Mn_at_stations):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and chargeaiBlity Mn
        they should be given on integration nodes.
        """
        sigma_oo = sigma_0 + Mn
        sigma_oo_faces = sigma_0_faces + Mn_faces
        sigma_oo_at_stations = sigma_0_at_stations + Mn_at_stations

        self.logger.debug("sigma_0 = %s" % (str(sigma_0)))
        self.logger.debug("Mn = %s" % (str(Mn)))
        self.logger.debug("sigma_oo = %s" % (str(sigma_oo)))

        self.logger.info("sigma_0 = %s" % (str(sigma_0)))
        self.logger.debug("sigma_0_face = %s" % (str(sigma_0_faces)))
        self.logger.debug("sigma_0_stations = %s - %s " % (min(sigma_0_at_stations), max(sigma_0_at_stations)))

        secondary_potentials_DC = getSecondaryPotentials(self.forward_pde,
                                                         sigma = sigma_0,
                                                         sigma_at_faces= sigma_0_faces,
                                                         schedule = self.data,
                                                         sigma_at_station = sigma_0_at_stations,
                                                         source_potential = self.source_potential,
                                                         sigma_src = self.sigma_src,
                                                         mask_faces = self.mask_outer_faces)
        secondary_potentials_DC_at_stations = {iA : self.grabValuesAtStations(secondary_potentials_DC[iA])
                                               for iA in secondary_potentials_DC}
        if self.logger.isEnabledFor(logging.DEBUG) :
            for iA in self.source_potential:
                self.logger.debug("DC secondary potential %d :%s" % (iA, str(secondary_potentials_DC[iA])))
        self.logger.info("%s DC secondary potentials calculated." % len(secondary_potentials_DC))

        potentials_S = getSecondaryPotentials(self.forward_pde,
                                              sigma = sigma_oo,
                                              sigma_at_faces= sigma_oo_faces,
                                              schedule = self.data,
                                              sigma_at_station = sigma_oo_at_stations,
                                              source_potential = self.source_potential,
                                              sigma_src = self.sigma_src,
                                              mask_faces = self.mask_outer_faces)
        for iA in potentials_S:
            potentials_S[iA]-=secondary_potentials_DC[iA]
        potentials_S_at_stations = {iA : self.grabValuesAtStations(potentials_S[iA])
                                    for iA in potentials_S}
        if self.logger.isEnabledFor(logging.DEBUG) :
            for iA in self.source_potential:
                self.logger.debug("IP potential %d :%s" % (iA, str(potentials_S[iA])))
        self.logger.info("%s IP potentials calculated." % len(potentials_S))

        return secondary_potentials_DC, secondary_potentials_DC_at_stations, potentials_S, potentials_S_at_stations

    def getMisfit(self, sigma_0, Mn, secondary_potentials_DC, secondary_potentials_DC_at_stations, potentials_S, potentials_S_at_stations, *args):
        """
        return the misfit in potential and chargeaiBlity weighted by weighting_misfit_DC factor
        """

        misfit_IP, misfit_DC = 0., 0.
        for iA, iB in self.misfit_IP:
            misfit_IP += self.misfit_IP[(iA, iB)].getValue(
                potentials_S_at_stations[iA] - potentials_S_at_stations[iB])
            misfit_DC += self.misfit_DC[(iA, iB)].getValue(
                 secondary_potentials_DC_at_stations[iA] - secondary_potentials_DC_at_stations[iB]
               + self.source_potentials_at_stations[iA] - self.source_potentials_at_stations[iB] )
        # print(misfit_IP, misfit_DC)
        return misfit_IP * self.weightingMisfit_2, misfit_DC * self.weightingMisfit_0

    def getDMisfit(self, sigma_0, secondary_potentials_DC, secondary_potentials_DC_at_stations, potentials_S, potentials_S_at_stations, *args):
        """
        returns the derivative of the misfit function with respect to sigma_0 and with respect to Mn
        """
        #TODO
        DMisfitDsigma_0 = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))
        DMisfitDMn = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))

        sigma_oo = (sigma_0 + Mn)


        #-------------------
        self.forward_pde.setValue(A=sigma_0 * kronecker(self.forward_pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for iA, iB in self.misfit_IP:
            A = self.data.getKeyOfStationNumber(iA)
            B = self.data.getKeyOfStationNumber(iB)
            dmis_2 = self.misfit_IP[(iA, iB)].getDerivative(
                Potentials_2_at_Stations[iA] - Potentials_2_at_Stations[iB])
            dmis_0 = self.misfit_DC[(iA, iB)].getDerivative(
                Potentials_0_at_Stations[iA] - Potentials_0_at_Stations[iB])

            Vs_star_2 = 0
            Q = 0
            for i in range(len(self.misfit_DC[(iA, iB)])):
                iM = self.misfit_DC[(iA, iB)].iMs[i]
                iN = self.misfit_DC[(iA, iB)].iNs[i]
                M = self.data.getKeyOfStationNumber(iM)
                N = self.data.getKeyOfStationNumber(iN)
                Vs_0 = self.injection_potential[M] - self.injection_potential[N] + Potentials_0[iM] - \
                       Potentials_0[iN]
                Vs_2 = Potentials_2[iM] - Potentials_2[iN]
                Vs_oo = Vs_0 - Vs_2

                Vs_star_2 += dmis_2[i] * Vs_oo
                Q += dmis_0[i] * Vs_0
            Vs_star_2 *= self.weightingMisfit_2
            self.forward_pde.setValue(
                X=sigma_0 * self.weightingMisfit_0 * grad(Q) + Mn * grad(Vs_star_2))  # point sources would be faster!
            Vs_star_0 = self.forward_pde.getSolution()

            Vs_0 = self.injection_potential[A] - self.injection_potential[B] + Potentials_0[iA] - Potentials_0[iB]
            Vs_2 = Potentials_2[iA] - Potentials_2[iB]
            Vs_oo = Vs_0 - Vs_2

            DMisfitDsigma_0 -= inner(grad(Vs_star_2), grad(Vs_2)) + inner(grad(Vs_star_0), grad(Vs_0))
            DMisfitDMn += inner(grad(Vs_star_2), grad(Vs_oo))

        return DMisfitDsigma_0, DMisfitDMn


class IPInversion(IPMisfitCostFunction):
    """
    Base class to run a IP inversion for conductivity (sigma_0, DC) and normalized chargeabilty (Mn=sigma_oo-sigma_0)
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1.e-4, Mn_ref=1.e-5, sigma_background=1.e-4,
                 w1=1., theta=0, usel1=False, epsl1=1e-4,
                 mask_fixed_property = None,
                 weighting_misfit_DC = 0.5, pde_tol=1.e-8, stationsFMT="e%s", logclip=5, logger=None, EPSILON=1e-15, **kargs):
        """       
        :domain: pde domain
        :data: survey data, is `fingal.SurveyData`
        :mask_fixed_property: mask of region where Mn and sigma_0 are fixed (set to reference values).
            If not set the bottom of the domain is used.
        :sigma_background: background conductivity, used to calculate the injection potentials
        :sigma_0_ref: reference DC conductivity
        :Mn_ref: reference normalized chargeability (must be positive)
        :param w1: weighting H1 regularization  int grad(m)^2
        :param theta: weighting cross-gradient term
        :weighting_misfit_DC: weighting factor for ERT part cost function, use 0.5
        :adjVs_ootStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format Vs_ooed to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        :EPSILON: absolute tolerance for zero values.
        """
        super().__init__(domain=domain, data=data, weighting_misfit_DC=weighting_misfit_DC, sigma_background=sigma_background,
                         pde_tol=pde_tol, stationsFMT=stationsFMT, logger=logger, **kargs)
        self.logclip = logclip
        self.EPS=EPSILON
        if mask_fixed_property is None:
            z = self.domain.getX()[2]
            self.mask_fixed_property = whereZero(z - sup(z))
        else:
            self.mask_fixed_property = mask_fixed_property
        self.sigma_0_ref = sigma_0_ref
        self.Mn_ref = Mn_ref
        self.usel1 = usel1
        self.epsl1 = epsl1


        # regularization
        self.w1 = w1
        self.theta=theta
        if self.usel1:
            self.Hpde = [ setupERTPDE(self.domain), setupERTPDE(self.domain) ]
            for pde in self.Hpde:
                pde.setValue(q=self.mask_fixed_property)
                pde.getSolverOptions().setTolerance(min(sqrt(pde_tol), 1e-3))
        else:
            self.Hpde = setupERTPDE(self.domain)
            self.Hpde.setValue(A=self.w1 * kronecker(3), q=self.mask_fixed_property)
            self.Hpde.getSolverOptions().setTolerance(min(sqrt(pde_tol), 1e-3))

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

    def getArguments(self, m):
        """
        precalculated parameters:
        """
        # recover temperature (interpolated to elements)
        im = interpolate(m, Function(self.domain))
        iMn = self.getMn(im[1])
        isigma_0 = self.getSigma0(im[0])
        args2 = self.getIPModelAndResponse(sigma_0=isigma_0, Mn=iMn)
        return im, isigma_0, iMn, args2

    def getValue(self, m, im, isigma_0, iMn, args2):
        """
        return the value of the cost function
        """
        misfit_2, misfit_0 = self.getMisfit(isigma_0, iMn, *args2)

        gm=grad(m, where=im.getFunctionSpace())
        if self.usel1:
            R=  self.w1 * ( integrate(sqrt(length(gm[0]) ** 2 + self.epsl1**2) ) + integrate(sqrt(length(gm[1]) ** 2 + self.epsl1**2) ) )
        else:
            R = self.w1 / 2 * integrate(length(gm) ** 2)

        if self.theta>0:
            gm0 = gm[0]
            gm1 = gm[1]
            lgm0 = length(gm0)
            lgm1 = length(gm1)

            lgm0 += whereZero(lgm0, tol=0) * self.EPS
            lgm1 += whereZero(lgm1, tol=0) * self.EPS
            MX = integrate(self.theta * ((lgm0 * lgm1) ** 2 - inner(gm0, gm1) ** 2) * (1 / lgm0 ** 2 + 1 / lgm1 ** 2))/2.
        else:
            MX=0
        V = MX + R + misfit_2 + misfit_0
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'misfit ERT, IP; reg ; Xgrad, total \t=  {misfit_0:e}, {misfit_2:e};  {R:e}; {MX:e} = {V:e}')
            self.logger.debug(
                f'ratios ERT, IP; reg ; Xgrad [%] \t=  {misfit_0/V*100:g}, {misfit_2/V*100:g};  {R/V*100:g}; {MX/V*100:g}')
        return V

    def getGradient(self, m, im, isigma_0, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gm=grad(m, where=im.getFunctionSpace())
        X = self.w1 * gm
        if self.usel1:
            X[0] *= 1/sqrt(length(gm[0]) ** 2 + self.epsl1 ** 2)
            X[1] *= 1/sqrt(length(gm[1]) ** 2 + self.epsl1 ** 2)

        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)

        Dsigma_0Dm0 = self.getDsigma0Dm(isigma_0, im[0])
        DMnDm1 = self.getDMnDm(iMn, im[1])
        Y=Data(0., (2,), im.getFunctionSpace())
        Y[1] = DMisfitDMn * DMnDm1
        Y[0] = DMisfitDsigma_0 * Dsigma_0Dm0

        if self.theta>0:
            gm0 = gm[0]
            gm1 = gm[1]
            lgm0 = length(gm0)
            lgm1 = length(gm1)
            lgm0 += whereZero(lgm0, tol=0) * self.EPS
            lgm1 += whereZero(lgm1, tol=0) * self.EPS

            X01 = inner(gm0, gm1)

            f = X01 * (1 / lgm0 ** 2 + 1 / lgm1 ** 2)
            X[0, :] += self.theta * ((1 + X01 ** 2 / lgm0 ** 4) * gm0 - f * gm1)
            X[1, :] += self.theta * ((1 + X01 ** 2 / lgm1 ** 4) * gm1 - f * gm0)

        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, m, im, isigma_0, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        if initializeHessian and self.usel1:
            gm=grad(m)
            for k in range(2):
                L=sqrt(length(gm[k]) ** 2 + self.epsl1 ** 2)
                self.Hpde[k].setValue(A=self.w1 * ( 1/L**3 * kronecker(3) - 1/L * outer(gm[k], gm[k])) )

        dm=Data(0., (2,), Solution(self.domain))
        if self.usel1:
            for i in [0,1]:
                self.Hpde[i].setValue(X=r[1][i], Y=r[0][i])
                dm[i] = self.Hpde[i].getSolution()
                txt = str(dm[i])
                if getMPIRankWorld() == 0:
                    self.logger.debug(f"search direction component {i} = {txt}.")
        else:
            for i in [0,1]:
                self.Hpde.setValue(X=r[1][i], Y=r[0][i])
                dm[i] = self.Hpde.getSolution()

                txt = str(dm[i])
                if getMPIRankWorld() == 0:
                    self.logger.debug(f"search direction component {i} = {txt}.")
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


class ERTMisfitCostFunction(CostFunction):
    """
    this the data misfit costfunction for IP
    """

    def __init__(self, domain=None, data=None, sigma_src=1., pde_tol=1.e-8,
                 mask_outer_faces = None, stationsFMT="e%s", useLogMisfit=False, logger=None, **kargs):
        """
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`
        :param weighting_misfit_ERT: weighting factor for chargaeaiBlity part cost function
        :stationsFMT: format string to map station keys A to mesh tags stationsFMT%A or None
        :logclip: cliping for p to avoid overflow in conductivity calculation
        """
        if logger is None:
            self.logger = logging.getLogger('fingal')
        else:
            self.logger = logger

        super().__init__(**kargs)
        self.domain = domain
        self.datatol = 1e-30
        self.stationsFMT = stationsFMT
        self.mask_outer_faces = makeMaskForOuterSurface(domain, facemask=mask_outer_faces)
        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.forward_pde = setupERTPDE(domain)
        self.forward_pde.getSolverOptions().setTolerance(pde_tol)
        self.data = data
        self.useLogMisfit=useLogMisfit
        # extract values  by the Locator
        station_locations = []
        for k in data.getStationNumeration():
            station_locations.append(data.getStationLocationByKey(k))
        #self.__grab_values_at_stations = Locator(Solution(domain), station_locations) # Dirac?
        #TODO: is this really working?
        self.__grab_values_at_stations = Locator(DiracDeltaFunctions(domain), station_locations)  # Dirac?

        self.sigma_src = sigma_src # needs to be a constant.
        self.setSourcePotentials()  # S_s is in the notes

        # build the misfit data (indexed by source index):
        self.misfit_DC = {}  # potential increment to injection field to get DC potential
        nd0 = 0 # counting number of data
        for A, B in self.data.injectionIterator():
            obs = self.data.getObservations(A, B)
            iA = self.data.getStationNumber(A)
            iB = self.data.getStationNumber(B)
            iMs = [self.data.getStationNumber(M) for M, N in obs]
            iNs = [self.data.getStationNumber(N) for M, N in obs]

            data_DC = np.array([self.data.getResistenceData((A, B, M, N)) for M, N in obs])
            #TODO NO ERROR IN WEIGHTING USED
            if self.useLogMisfit:
                error_DC = max ([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs])
                self.misfit_0[(iA, iB)] =  DataMisfitLog(iMs=iMs, data=data_0, iNs=iNs, injections=(A,B), weightings=1. / error_DC**2)
            else:
                error_DC = max ([self.data.getResistenceError((A, B, M, N)) for M, N in obs])
                self.misfit_DC[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_DC, iNs=iNs, injections=(A, B), weightings =1./error_DC**2/data_DC.size)

            nd0+= len(self.misfit_DC[(iA, iB)])
        self.logger.info("%d DC data records detected." % (nd0))

        if nd0 > 0:
            for iA, iB in self.misfit_DC:
                self.misfit_DC[(iA, iB)].rescaleWeight(1. / nd0)
    def grabValuesAtStations(self, m):
        """
        returns the values of m at the stations. The order is given by data.getStationNumeration()
        """
        return np.array(self.__grab_values_at_stations(m))
    def getObservationsFromSourcePotential(self, A  , B ):
        """
        calculate observations MN for injection AB for the source potentials.
        ordering follows data.getObservations(A, B)
        """
        obs = self.data.getObservations(A, B)
        iA = self.data.getStationNumber(A)
        iB = self.data.getStationNumber(B)
        injectionsAB= self.source_potentials_at_stations[iA] - self.source_potentials_at_stations[iB]
        source_obs=np.array([injectionsAB[self.data.getStationNumber(M)] - injectionsAB[self.data.getStationNumber(N)] for M, N in obs])
        return source_obs
    def setNewSigmaSrc(self, new_sigma_src):
        """
        updated sigma_background
        """
        factor= self.sigma_src / new_sigma_src
        for iA in self.source_potential:
            self.source_potential[iA]*=factor
            self.source_potentials_at_stations[iA]*=factor
        self.sigma_src=new_sigma_src
    def setSourcePotentials(self):
        """
        return the primary electric potential for all injections (A,B) Vs_ooing conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_Field
        """
        self.logger.debug("source conductivity sigma_src =  %s" % (str(self.sigma_src)))
        self.source_potential = getSourcePotentials(self.domain, self.sigma_src, self.data, mask_outer_faces =self.mask_outer_faces, stationsFMT=self.stationsFMT)
        self.source_potentials_at_stations = { iA : self.grabValuesAtStations(self.source_potential[iA])
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
        if self.useLogMisfit:
            for iA, iB in self.misfit_DC:
                u = self.source_potentials_at_stations[iA] - self.source_potentials_at_stations[iB]
                d = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * (self.misfit_DC[(iA, iB)].data - log(abs(d)) ) )
                b+= self.misfit_DC[(iA, iB)].data.size
            if b>0:
                beta=exp(a/b)
        else:
            for iA, iB in self.misfit_DC:
                u = self.source_potentials_at_stations[iA] - self.source_potentials_at_stations[iB]
                d = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * self.misfit_DC[(iA, iB)].data * d )
                b+= sum(  self.misfit_DC[(iA, iB)].weightings * d**2 )
            if b > 0:
                beta = a/b

        assert beta > 0, "conductivity correction factor must be positive."
        return beta * self.sigma_src

    def getERTModelAndResponse(self, sigma_0, sigma_0_face, sigma_0_at_stations):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and sigma_0_face
        they should be given on integration nodes.
        """
        self.logger.info("sigma_0 = %s" % (str(sigma_0)))
        self.logger.debug("sigma_0_face = %s" % (str(sigma_0_face)))
        self.logger.debug("sigma_0_stations = %s - %s " % (min(sigma_0_at_stations), max(sigma_0_at_stations)))

        secondary_potentials_DC = getSecondaryPotentials(self.forward_pde,
                                                         sigma = sigma_0,
                                                         sigma_at_faces= sigma_0_face,
                                                         schedule = self.data,
                                                         sigma_at_station = sigma_0_at_stations,
                                                         source_potential = self.source_potential,
                                                         sigma_src = self.sigma_src,
                                                         mask_faces = self.mask_outer_faces)
        secondary_potentials_DC_at_stations = {iA : self.grabValuesAtStations(secondary_potentials_DC[iA])
                                               for iA in secondary_potentials_DC}
        if self.logger.isEnabledFor(logging.DEBUG) :
            for iA in self.source_potential:
                self.logger.debug("DC secondary potential %d :%s" % (iA, str(secondary_potentials_DC[iA])))
        self.logger.info("%s DC secondary potentials calculated." % len(secondary_potentials_DC))
        return secondary_potentials_DC, secondary_potentials_DC_at_stations

    def getMisfit(self, sigma_0, sigma_0_face, sigma_0_stations, secondary_potentials_0, secondary_potentials_0_at_stations, *args):
        """
        return the misfit in potential and chargeaiBlity weighted by weighting_misfit_DC factor
        """
        potentials_DC_at_stations = {}
        for iA in secondary_potentials_0_at_stations:
            potentials_DC_at_stations[iA] = self.source_potentials_at_stations[iA] + secondary_potentials_0_at_stations[iA]
        misfit_0 = 0.
        for iA, iB in self.misfit_DC:
            misfit_0 += self.misfit_DC[(iA, iB)].getValue(
                potentials_DC_at_stations[iA] - potentials_DC_at_stations[iB])
        return misfit_0
    def getDMisfit(self, sigma_0, sigma_0_face, sigma_0_stations, secondary_potentials_DC, secondary_potentials_0_at_stations, *args):
        """
        returns the derivative of the misfit function with respect to sigma_0 with respect to Mn
        """
        SOURCES=np.zeros( (self.data.getNumStations(), self.data.getNumStations()), float)
        potentials_0_at_stations = {}
        for iA in secondary_potentials_0_at_stations:
            potentials_0_at_stations[iA] = self.source_potentials_at_stations[iA] + secondary_potentials_0_at_stations[iA]

        for iA, iB in self.misfit_DC:
            dmis_0 = self.misfit_DC[(iA, iB)].getDerivative(
                potentials_0_at_stations[iA] - potentials_0_at_stations[iB])
            for i in range(len(self.misfit_DC[(iA, iB)])):
                iM = self.misfit_DC[(iA, iB)].iMs[i]
                iN = self.misfit_DC[(iA, iB)].iNs[i]
                #ABMN
                SOURCES[iA, iM] += dmis_0[i]
                SOURCES[iB, iM] -= dmis_0[i]
                SOURCES[iA, iN] -= dmis_0[i]
                SOURCES[iB, iN] += dmis_0[i]

        DMisfitDsigma_0 = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))
        DMisfitDsigma_0_face = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('y'))
        self.forward_pde.setValue(A=sigma_0 * kronecker(self.forward_pde.getDim()), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
        n = self.domain.getNormal()
        x = FunctionOnBoundary(self.domain).getX()

        for iA in secondary_potentials_DC.keys():
                s = Scalar(0., DiracDeltaFunctions(self.forward_pde.getDomain()))
                for M in self.data.getStationNumeration():
                    iM = self.data.getStationNumber(M)
                    if self.stationsFMT is None:
                        s.setTaggedValue(M, SOURCES[iA, iM])
                    else:
                        s.setTaggedValue(self.stationsFMT % M,  SOURCES[iA, iM])
                self.forward_pde.setValue(y_dirac=s)
                r = x - self.data.getStationLocationByNumber(iA)
                fA = inner(r, n) / length(r) ** 2 * self.mask_outer_faces
                self.forward_pde.setValue(d=sigma_0_face * fA )
                VA_star=self.forward_pde.getSolution()
                self.logger.debug("DC adjoint potential %d :%s" % (iA, str(VA_star)))
                VA= self.source_potential[iA] + secondary_potentials_DC[iA]
                DMisfitDsigma_0 -= inner(grad(VA_star), grad(VA))
                DMisfitDsigma_0_face -= ( fA * VA_star ) * VA

        self.logger.info("%s adjoint potentials calculated." % len(secondary_potentials_DC.keys()))
        return DMisfitDsigma_0, DMisfitDsigma_0_face


class ERTInversion(ERTMisfitCostFunction):
    """
    Base class to run a IP inversion for conductivity (sigma_0, DC)
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1., useL1Norm=False, epsilonL1Norm=1e-4,
                 mask_fixed_property = None, mask_outer_faces = None,
                 pde_tol=1.e-8, reg_tol=None, stationsFMT="e%s", logclip=5,
                 useLogMisfit=False, logger=None, EPSILON=1e-15, **kargs):
        """
        :param domain: pde domain
        :param data: survey data, is `fingal.SurveyData`
        :param sigma_0_ref: sigma = sigma_0_ref * exp(m)
        :param sigma_src: conductivity used to calculate the injection/source potentials (need to be constant)
        :param w1: weighting H1 regularization  int grad(m)^2
        :param useL1Norm: use L1+eps regularization
        :param epsilonL1Norm: cut-off for L1+eps regularization
        :param mask_fixed_property: mask of region where sigma_0 is fixed to sigma_0_ref.
            If not set the bottom of the domain is used.
        :param mask_outer_faces: mask of the faces where potential field 'radiation' is set.
        :param pde_tol: tolerance for the PDE solver
        :param reg_tol: tolerance for the PDE solver in the regularization
        :param stationsFMT: format to map station keys k to mesh tags stationsFMT%k or None
        :param logclip: cliping for p to avoid overflow in conductivity calculation
        :param EPSILON: absolute tolerance for zero values.
        """
        if sigma_src == None:
            sigma_src = sigma_0_ref
        super().__init__(domain=domain, data=data, sigma_src=sigma_src, pde_tol=pde_tol,
                         mask_outer_faces = makeMaskForOuterSurface(domain,mask_outer_faces ),
                         stationsFMT=stationsFMT, logger=logger, useLogMisfit=useLogMisfit, **kargs)
        self.logclip = logclip
        self.EPS=EPSILON
        # its is assumed that sigma_ref on the gaces is not updated!!!
        if mask_fixed_property is None:
            z = self.forward_pde.getDomain().getX()[2]
            self.mask_fixed_property  =  whereZero(z - inf(z))
        else:
            self.mask_fixed_property = wherePositive( mask_fixed_property)
        self.sigma_0_ref = sigma_0_ref
        self.useL1Norm = useL1Norm
        self.epsilonL1Norm = epsilonL1Norm
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))


        # regularization
        self.w1 = w1
        self.factor2D=1
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.setValue(q=self.mask_fixed_property)
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.HpdeUpdateCount=2
        K=kronecker(3)
        K[1,1]*=self.factor2D
        self.Hpde.setValue(A=self.w1 * K)

        #  reference conductivity:
        self.setSigma0Ref(sigma_0_ref)


    def setSigma0Ref(self, sigma_0_ref):
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

    def getArguments(self, m):
        """
        precalculated parameters:
        """
        # recover temperature (interpolated to elements)
        im = interpolate(m, Function(self.domain))
        im_face =interpolate(m, FunctionOnBoundary(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % ( str(im)))

        isigma_0 = self.getSigma0(im)
        isigma_0_face = self.getSigma0(im_face)
        isigma_0_stations =  self.getSigma0(im_stations)
        args2 = self.getERTModelAndResponse(isigma_0, isigma_0_face, isigma_0_stations)
        return im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2

    def getValue(self, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)

        gm=grad(m, where=im.getFunctionSpace())
        lgm2=gm[0]**2+self.factor2D*gm[1]**2+gm[2]**2
        if self.useL1Norm:
            R=  self.w1 * integrate(sqrt(lgm2 + self.epsilonL1Norm ** 2))
        else:
            R = self.w1 / 2 * integrate(lgm2)

        V = R + misfit_0
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'misfit ERT; reg ; total \t=  {misfit_0:e}, ;{R:e}; = {V:e}')
            self.logger.debug(
                f'ratios ERT; reg  [%] \t=  {misfit_0/V*100:g}; {R/V*100:g}')
        return V

    def getGradient(self, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gm=grad(m, where=im.getFunctionSpace())
        X =gm
        X[1]*=self.factor2D
        if self.useL1Norm:
            lgm2 = gm[0] ** 2 + self.factor2D * gm[1] ** 2 + gm[2] ** 2
            X *= 1/sqrt(lgm2 + self.epsilonL1Norm ** 2)
        X*=self.w1
        DMisfitDsigma_0,  DMisfitDsigma_0_face = self.getDMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)

        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Dsigma_0Dm_face = self.getDsigma0Dm(isigma_0_face, im_face)
        Y = DMisfitDsigma_0 * Dsigma_0Dm
        y = DMisfitDsigma_0_face * Dsigma_0Dm_face
        return ArithmeticTuple(Y, X, y)

    def getInverseHessianApproximation(self, r, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        clip_property_function=10
        if self.useL1Norm:
            if self.HpdeUpdateCount > 1:
                gm=grad(m)
                lgm2 = gm[0] ** 2 + self.factor2D * gm[1] ** 2 + gm[2] ** 2
                L = sqrt(lgm2 + self.epsilonL1Norm ** 2)
                self.Hpde.setValue(A=self.w1 * ( 1/L  * kronecker(3) - 1/L**3 * outer(gm, gm)) )
                self.HpdeUpdateCount+=1
                print("TO DO")

        self.Hpde.setValue(X=r[1], Y=r[0], y=r[2])
        dm = self.Hpde.getSolution()
        #dm.copyWithMask(clip_property_function+0*m, wherePositive(m+dm-clip_property_function))
        #dm.copyWithMask(-clip_property_function+0*m, wherePositive(-clip_property_function-(m+dm)))
        self.logger.debug(f"search direction component = {str(dm)}.")
        return dm
    def getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(r[0] * m + inner(r[1], grad(m)) ) + integrate(r[2] * m )

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)