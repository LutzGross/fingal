"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException
from .tools import setupERTPDE, getSourcePotentials, makeMaskForOuterSurface, getSecondaryPotentials, DataMisfitQuad, DataMisfitLog
import logging
import numpy as np
from esys.escript.pdetools import Locator,  ArithmeticTuple


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