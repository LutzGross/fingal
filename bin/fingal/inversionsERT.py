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
                #error_DC = max ([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs])
                error_DC = np.array([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs])
                self.misfit_DC[(iA, iB)] =  DataMisfitLog(iMs=iMs, data=data_DC, iNs=iNs, injections=(A,B), weightings=1. / error_DC**2/data_DC.size)
            else:
                #error_DC = max ([self.data.getResistenceError((A, B, M, N)) for M, N in obs])
                error_DC = np.array([self.data.getResistenceError((A, B, M, N)) for M, N in obs])
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
        self.source_potential = getSourcePotentials(self.domain, self.sigma_src, self.data, mask_outer_faces =self.mask_outer_faces, stationsFMT=self.stationsFMT, logger=self.logger)
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
                du = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * (self.misfit_DC[(iA, iB)].data - log(abs(du)) ))
                b+= sum(  self.misfit_DC[(iA, iB)].weightings )
            if b>0:
                beta=exp(a/b)
        else:
            for iA, iB in self.misfit_DC:
                u = self.source_potentials_at_stations[iA] - self.source_potentials_at_stations[iB]
                du = self.misfit_DC[(iA, iB)].getDifference(u)
                a+= sum(  self.misfit_DC[(iA, iB)].weightings * self.misfit_DC[(iA, iB)].data * du )
                b+= sum(  self.misfit_DC[(iA, iB)].weightings * du**2 )
            if b > 0:
                beta = a/b

        assert beta > 0, "conductivity correction factor must be positive."
        return self.sigma_src/beta

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
                                                         mask_faces = self.mask_outer_faces, logger=self.logger)
        secondary_potentials_DC_at_stations = {iA : self.grabValuesAtStations(secondary_potentials_DC[iA])
                                               for iA in secondary_potentials_DC}
        #if self.logger.isEnabledFor(logging.DEBUG) :
        #    for iA in self.source_potential:
        #        self.logger.debug("DC secondary potential %d :%s" % (iA, str(secondary_potentials_DC[iA])))
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
                #self.logger.debug("DC adjoint potential %d :%s" % (iA, str(VA_star)))
                VA= self.source_potential[iA] + secondary_potentials_DC[iA]
                DMisfitDsigma_0 -= inner(grad(VA_star), grad(VA))
                DMisfitDsigma_0_face -= ( fA * VA_star ) * VA

        self.logger.debug("%s adjoint potentials calculated." % len(secondary_potentials_DC.keys()))
        return DMisfitDsigma_0, DMisfitDsigma_0_face


class ERTInversionH1(ERTMisfitCostFunction):
    """
    ERT inversion with H1 regularization
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
            x = self.forward_pde.getDomain().getX()
            qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
            qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
            qz = whereZero(x[2] - inf(x[2]))
            self.mask_fixed_property  =  0 * qx + 0* qy + qz
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
        self.updateSigma0Ref(sigma_0_ref)


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
    def extractPropertyFunction(self, m):
        return m
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

class ERTInversionPseudoGauss(ERTMisfitCostFunction):
    """
    ERT inversion with Gaussian regulization with decoupled preconditioning
    """
    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1., length_scale =1.,
                mask_outer_faces = None,
                 pde_tol=1.e-8, reg_tol=None, stationsFMT="e%s", logclip=5,
                 useLogMisfit=False, logger=None, EPSILON=1e-15, **kargs):
        """
        :param domain: pde domain
        :param data: survey data, is `fingal.SurveyData`
        :param sigma_0_ref: sigma = sigma_0_ref * exp(m)
        :param sigma_src: conductivity used to calculate the injection/source potentials (need to be constant)
        :param w1: weighting H2 regularization int (m+a^2*del(m)) ^2
        :param length_scale: length scale factor a
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
        self.sigma_0_ref = sigma_0_ref
        self.length_scale = length_scale
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))
        self.logger.info("length scale = %s" % (str(self.length_scale)))

        # regularization
        self.w1 = w1
        self.Hpde = setupERTPDE(self.domain)
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.HpdeUpdateCount=1
        self.Hpde.setValue(A=self.w1 * self.length_scale**2 * kronecker(3), D=self.w1 )

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

    def extractPropertyFunction(self, M):
        return M[0]
    def getArguments(self, M):
        """
        precalculated parameters:
        """
        m=M[0]
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

    def getValue(self, M, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)
        a=self.length_scale

        gradm =grad(M[0], where=im.getFunctionSpace())
        divM = div(M[1:], where=im.getFunctionSpace())
        curlM = curl(M[1:], where = im.getFunctionSpace())

        reg = self.w1 * integrate(length(im-a*divM)**2)
        gM =  self.w1 * integrate(length(M[1:]-a*gradm)**2)
        cM =  self.w1 * integrate(length(a * curlM)**2)

        R = (reg + gM + cM )/2
        V = R + misfit_0
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'misfit ERT; reg; div; curl; total \t=  {misfit_0:e} ;{reg:e}; {gM:e}; {cM:e};= {V:e}')
            self.logger.debug(
                f'ratios ERT; reg  [%] \t=  {misfit_0/V*100:g}; {R/V*100:g}')
        return V

    def getGradient(self, M, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        a=self.length_scale

        gradm =grad(M[0], where=im.getFunctionSpace())
        divM = div(M[1:], where=im.getFunctionSpace())
        curlM = curl(M[1:], where = im.getFunctionSpace())

        X = Data(0, (4,3), Function(M.getDomain()))
        Y = Data(0, (4,),  Function(M.getDomain()))
        y = Data(0, (4,),  FunctionOnBoundary(M.getDomain()))

        Y[1:]+= M[1:]- a * gradm
        X[0,:] = - a * (M[1:]- a * gradm)

        d=M[0]-a*divM
        Y[0]+=d
        for i in range(0, self.domain.getDim()):
                X[1+i,i]= - a * d

        X[3,1] = a**2 * curlM[0]
        X[2,2] = - X[3,1]
        X[1,2] = a**2 * curlM[1]
        X[3,0] = - X[1,2]
        X[2,0] =  a**2 * curlM[2]
        X[1,1] = - X[2,0]

        X*=self.w1
        Y*=self.w1
        #
        DMisfitDsigma_0,  DMisfitDsigma_0_face = self.getDMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)
        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Dsigma_0Dm_face = self.getDsigma0Dm(isigma_0_face, im_face)
        Y[0] += DMisfitDsigma_0 * Dsigma_0Dm
        y[0] += DMisfitDsigma_0_face * Dsigma_0Dm_face

        return ArithmeticTuple(Y, X, y)

    def getInverseHessianApproximation(self, r, M, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        dM=Data(0., M.getShape(), M.getFunctionSpace() )
        for i in range(M.getDomain().getDim()+1):
            self.Hpde.setValue(X=r[1][i,:], Y=r[0][i], y=r[2][i])
            dM[i] = self.Hpde.getSolution()
            self.logger.debug(f"search direction component {i} = {str(dM[i])}.")
        return dM
    def getDualProduct(self, M, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(inner(r[0],M) + inner(r[1], grad(M))) + integrate(inner(r[2], M) )

    def getNorm(self, M):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(M)

class ERTInversionH2(ERTMisfitCostFunction):
    """
    ERT inversion with H2 regularization
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1.,
                 mask_outer_faces = None,
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
        self.sigma_0_ref = sigma_0_ref
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))
        # masks for surfaces:
        x=self.forward_pde.getDomain().getX()
        qx = whereZero(x[0] - inf(x[0])) + whereZero(x[0] - sup(x[0]))
        qy = whereZero(x[1] - inf(x[1])) + whereZero(x[1] - sup(x[1]))
        qz = whereZero(x[2] - inf(x[2]))


        # regularization
        self.mpde = setupERTPDE(self.domain)
        self.mpde.getSolverOptions().setTolerance(pde_tol)
        self.mpde.setValue(A=kronecker(3), q= 0*qx + 0*qy + qz )

        # regularization
        self.w1 = w1
        self.Hpde = setupERTPDE(self.domain)
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.HpdeUpdateCount=2
        self.Hpde.setValue(A=self.w1 * kronecker(3))
        self.Hpde_qs = [ qy + qz, qx + qz,  qx + qy ]
        #  reference conductivity:
        self.updateSigma0Ref(sigma_0_ref)


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

    def extractPropertyFunction(self, M):
        self.mpde.setValue(X=interpolate(M, Function(M.getDomain())), Y=Data(), y=Data())
        m=self.mpde.getSolution()
        return m
    def getArguments(self, M):
        """
        precalculated parameters:
        """
        m=self.extractPropertyFunction(M)

        im = interpolate(m, Function(self.domain))
        im_face =interpolate(m, FunctionOnBoundary(self.domain))
        im_stations = self.grabValuesAtStations(m)
        self.logger.debug("m = %s" % ( str(im)))

        isigma_0 = self.getSigma0(im)
        isigma_0_face = self.getSigma0(im_face)
        isigma_0_stations =  self.getSigma0(im_stations)
        args2 = self.getERTModelAndResponse(isigma_0, isigma_0_face, isigma_0_stations)
        return m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2

    def getValue(self, M, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        return the value of the cost function
        """
        misfit_0 = self.getMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)

        gM=grad(M, where=im.getFunctionSpace())
        R = self.w1 / 2 * integrate(length(gM)**2)

        V = R + misfit_0
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'misfit ERT; reg ; total \t=  {misfit_0:e}, ;{R:e}; = {V:e}')
            self.logger.debug(
                f'ratios ERT; reg  [%] \t=  {misfit_0/V*100:g}; {R/V*100:g}')
        return V

    def getGradient(self, M, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        gM=grad(M, where=im.getFunctionSpace())
        X = gM * self.w1
        DMisfitDsigma_0,  DMisfitDsigma_0_face = self.getDMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)

        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Dsigma_0Dm_face = self.getDsigma0Dm(isigma_0_face, im_face)
        self.mpde.setValue(X=Data(), Y = DMisfitDsigma_0 * Dsigma_0Dm, y = DMisfitDsigma_0_face * Dsigma_0Dm_face)
        Y=grad(self.mpde.getSolution(), where=X.getFunctionSpace())
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, M, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        dM = Vector(0., M.getFunctionSpace())
        for k in range(3):
            if k==10:
                dM[k]=0
            else:
                self.Hpde.setValue(q= self.Hpde_qs[k], X=r[1][k,:], Y=r[0][k])
                dM[k] = self.Hpde.getSolution()
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

class ERTInversionGauss(ERTMisfitCostFunction):
    """
    ERT inversion with H2 regularization
    """

    def __init__(self, domain=None, data=None,
                 sigma_0_ref=1e-4, sigma_src=None, w1=1.,  length_scale =1.,
                 mask_outer_faces = None,
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
        self.sigma_0_ref = sigma_0_ref
        self. length_scale = length_scale
        self.logger.info("Reference conductivity sigma_0_ref = %s" % (str(self.sigma_0_ref)))
        self.logger.info("length scale = %s" % (str(self.length_scale)))

        # regularization
        self.mpde = setupERTPDE(self.domain)
        self.mpde.getSolverOptions().setTolerance(pde_tol)
        self.mpde.setValue(A=self.length_scale**2 * kronecker(3), D=1)

        # regularization
        self.w1 = w1
        self.Hpde = setupERTPDE(self.domain)
        if not reg_tol:
            reg_tol=min(sqrt(pde_tol), 1e-3)
        reg_tol = pde_tol
        self.logger.debug(f'Tolerance for solving regularization PDE is set to {reg_tol}')
        self.Hpde.getSolverOptions().setTolerance(reg_tol)
        self.HpdeUpdateCount=2
        self.Hpde.setValue(D=self.w1 * 1.)
        #  reference conductivity:
        self.updateSigma0Ref(sigma_0_ref)


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

    def extractPropertyFunction(self, m):
        self.mpde.setValue(Y=interpolate(m, Function(m.getDomain())), X=Data(), y=Data())
        m=self.mpde.getSolution()
        return m
    def getArguments(self,m):
        """
        precalculated parameters:
        """
        m=self.extractPropertyFunction(m)

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

        R = self.w1 / 2 * integrate(im**2)

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
        Y = im * self.w1
        DMisfitDsigma_0,  DMisfitDsigma_0_face = self.getDMisfit(isigma_0, isigma_0_face, isigma_0_stations, *args2)

        Dsigma_0Dm = self.getDsigma0Dm(isigma_0, im)
        Dsigma_0Dm_face = self.getDsigma0Dm(isigma_0_face, im_face)
        self.mpde.setValue(X=Data(), Y = DMisfitDsigma_0 * Dsigma_0Dm, y = DMisfitDsigma_0_face * Dsigma_0Dm_face)
        Y+=self.mpde.getSolution()
        return Y

    def getInverseHessianApproximation(self, r, m, im, im_face, isigma_0, isigma_0_face, isigma_0_stations, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        self.Hpde.setValue(Y=r)
        dm = self.Hpde.getSolution()
        self.logger.debug(f"search direction = {str(dm)}.")
        return dm

    def getDualProduct(self, m, r):
        """
        dual product of gradient `r` with increment `m`. Overwrites `getDualProduct` of `MeteredCostFunction`
        """
        return integrate(r * m)
    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)
