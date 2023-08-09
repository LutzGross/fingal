"""
cost functions for ERT/IP inversions

by l.gross@uq.edu.au, 2021, 2028
"""

from esys.escript import *
from esys.escript.minimizer import CostFunction, MinimizerException
from .tools import setupERTPDE, getInjectionPotentials


class DataMisfitQuad(object):
    """
    A mechanism to quadratic handel misfits 
    
    sum weighting[i]*(1-observations[i]/data[i])**2  
    
    where observations are collected from a value vector iMs or difference using iMs and iNs
    """
    def __init__(self, iMs, data, iNs=None, weightings=None, data_tol=1e-8):
        self.Ms = iMs
        self.Ns = iNs
        self.data_tol = data_tol * np.ones(data.shape)
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
        return len(self.iMs)

    def getValue(self, u):
        if self.iNs:
            diff = u[self.iMs] - u[self.iNs]
        else:
            diff = u[self.iMs]

        R = diff / self.data
        res = 1 - R
        dd = res ** 2 * self.weightings
        return 0.5 * sum(dd)

    def getWeightedDifference(self, u):
        if self.iNs:
            diff = u[self.iMs] - u[self.iNs]
        else:
            diff = u[self.iMs]
        R = diff / self.data
        res = 1 - R
        dd = -res / self.data * self.weightings
        return dd

    def rescaleWeight(self, weight=1.):
        self.weightings *= weight

class DataMisfitLog2(object):
    """
    A mechanism to handel a logarithmic misfit 

    sum_i weighting[i]*(log((observations[i]/data[i] - 2) * observations[i]/data[i] + 2))**2  

    where observations are collected from a value vector iMs or difference using iMs and iNs
    """
    def __init__(self, iMs, data, iNs=None, weightings=None, data_tol=1e-8):
        self.iMs = iMs
        self.iNs = iNs
        self.data_tol = data_tol * np.ones(data.shape)
        # assert not ( weightings is None or quadWeight is None )

        if isinstance(data, np.ndarray):
            self.data = data
        else:
            self.data = np.array(data)

        if isinstance(weightings, np.ndarray):
            self.weightings = weightings
        elif weightings is None:
            self.weightings = np.zeros(self.data.shape)
        elif isinstance(weightings, float) or isinstance(weightings, int):
            self.weightings = np.zeros(self.data.shape) + weightings
        else:
            self.weightings = np.array(weightings)

    def __len__(self):
        return len(self.iMs)

    def getValue(self, u):
        if self.iNs:
            diff = u[self.iMs] - u[self.iNs]
        else:
            diff = u[self.iMs]

        R = diff / self.data
        res = log((R - 2) * R + 2)
        dd = res ** 2 * self.weightings
        # print("DATA LOG: data=", self.data, "result = ", diff, "ratio = ", R, "misfit =",res, "->", dd, 'sum', sum(dd))
        return 0.5 * sum(dd)

    def getWeightedDifference(self, u):
        if self.iNs:
            diff = u[self.iMs] - u[self.iNs]
        else:
            diff = u[self.iMs]
        R = diff / self.data
        res = log((R - 2) * R + 2)
        dd = 2 * res / ((R - 2) * R + 2) * (R - 1) / self.data * self.weightings
        return dd

    def rescaleWeight(self, weight=1.):
        self.weightings *= weight




class IPMisfitCostFunction(CostFunction):
    """
    this the data misfit costfunction for IP
    """

    def __init__(self, domain=None, data=None, weighting_misfit_0=0.5, sigma_background=1.e-05, pde_tol=1.e-8,
                 stationsFMT="e%s", logger=None, **kargs):
        """       
        :param domain: pde domain
        :param data: data, is `fingal.SurveyData`
        :param weighting_misfit_0: weighting factor for chargaeaiBlity part cost function
        :adjVs_ootStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format Vs_ooed to map station keys k to mesh tags stationsFMT%k or None
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
        self.setMisfitWeighting(weighting_misfit_0)

        # setup PDE for forward models (potentials are fixed on all faces except the surface)
        self.forward_pde = setupERTPDE(domain)
        self.forward_pde.getSolverOptions().setTolerance(pde_tol)
        # self.forward_pde.getSolverOptions().setTrilinosParameter("reVs_ooe: type","full")
        # self.forward_pde=setupERTPDE(domain)
        x = self.forward_pde.getdomain().getX()[0]
        y = self.forward_pde.getdomain().getX()[1]
        z = self.forward_pde.getdomain().getX()[2]
        q = whereZero(x - inf(x)) + whereZero(x - sup(x)) + whereZero(y - inf(y)) + whereZero(y - sup(y)) + whereZero(
            z - inf(z))

        self.forward_pde.setValue(q=q)
        self.data = data

        # set up stations in a specific order so one can work with arrays D
        # extracted by the Locator
        # D[i] = value at station M for i=getStationNumber(M) or M=getKeyOfStationNumber(i)
        station_locations = []
        for s in data.getStationNumeration():
            station_locations.append(data.getStationLocation(s))
        self.grabValues = Locator(Solution(domain), station_locations)
        # self.grabValues=Locator(DiracDeltaFunctions(domain), station_locations)

        self.sigma_background = interpolate(sigma_background, self.forward_pde.getFunctionSpaceForCoefficient('A'))
        self.setInjectionPotentials()  # S_s is in the notes

        # build the misfit data (indexed by source index):
        self.misfit_2 = {}  # secondary potential
        self.misfit_0 = {}  # potential increment to injection field to get DC potential
        nd0 = 0 # counting number of data
        nd2 = 0
        for A, B in self.data.injectionIterator():
            obs = self.data.getObservations(A, B)
            iA = self.data.getStationNumber(A)
            iB = self.data.getStationNumber(B)
            #print(A,B,"->",iA, iB)
            # ====
            # print(" --- ", A,B, " --- ")
            # for M, N in obs:
            #  print(M, N, self.data.getResistenceError((A,B,M,N))/self.data.getResistenceData((A,B,M,N)), "- ", 
            #        self.data.getSecondaryResistenceError((A,B,M,N))/self.data.getSecondaryResistenceData((A,B,M,N))  )
            # ====
            iMs = [self.data.getStationNumber(M) for M, N in obs]
            iNs = [self.data.getStationNumber(N) for M, N in obs]
            injections= self.injection_potential_at_stations[A] - self.injection_potential_at_stations[B]

            data_0 = np.array([self.data.getResistenceData((A, B, M, N)) for M, N in obs])  # u_oo
            rel_error_0 = np.array([self.data.getResistenceRelError((A, B, M, N)) for M, N in obs])
            # print('Data RES =', data_0, 1/rel_tol_0**2)

            # injection potentials are removed from the data:
            data_0 = data_0 - [injections[self.data.getStationNumber(M)] - injections[self.data.getStationNumber(N)] for
                               M, N in obs]
            # self.misfit_0[(iA,iB)]=DataMisfitLog2( Ms=Ms, data=data_0, Ns=Ns, weightings=1/rel_tol_0**2) #  test against Ws
            self.misfit_0[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_0, iNs=iNs,
                                                     weightings=1. / rel_error_0 ** 2)  # test against Ws
            nd0 += len(self.misfit_0[(iA, iB)])

            data_2 = []
            data_2_w = []
            for M, N in obs:
                v = self.data.getSecondaryResistenceData((A, B, M, N))  # u_0-u_oo
                e = getSecondaryResistenceError((A, B, M, N))
                if self.data.isUndefined(v):
                    data_2.append(-1000.)
                    data_2_w.append(0.)
                else:

                    data_2.append(v)
                    data_2_w.append( (v/e) ** 2)
                    nd2 += 1
            data_2, data_2_w = np.array(data_2), np.array(data_2_w)
            # self.misfit_2[(iA,iB)]=DataMisfitLog2( Ms=Ms, data=data_2, Ns=Ns, weightings=data_2_w) # test against Vs
            self.misfit_2[(iA, iB)] = DataMisfitQuad(iMs=iMs, data=data_2, iNs=iNs, weightings=data_2_w)  # test against Vs

        if getMPIRankWorld() == 0:
            self.logger.info("%d/%d DC/IP data records detected." % (nd2, nd0))

        nd = nd0 + nd2
        if nd > 0:
            for iA, iB in self.misfit_2:
                self.misfit_0[(iA, iB)].rescaleWeight(1. / nd)
                self.misfit_2[(iA, iB)].rescaleWeight(1. / nd)

    def setMisfitWeighting(self, weighting_misfit_0=0):

        self.weightingMisfit_2 = max(0, 1 - weighting_misfit_0)
        self.weightingMisfit_0 = max(0, weighting_misfit_0)

        # self.weightingMisfit_2=0
        # self.weightingMisfit_0=1

        # assert self.weightingMisfit_2>0 or self.weightingMisfit_0>0

    def setInjectionPotentials(self):
        """
        return the primary electric potential for all injections (A,B) Vs_ooing conductivity sigma
        :sigma: (primary) conductivity distribution
        :return: dictonary of injections (A,B)->primary_Field
        """
        txt = str(self.sigma_background)
        if getMPIRankWorld() == 0:
            self.logger.debug("reference conductivity =  %s" % (txt))
        self.injection_potential = getInjectionPotentials(self.domain, self.sigma_background, self.data, stationsFMT=self.stationsFMT)
        self.injection_potential_at_stations = {}

        for A in self.injection_potential:
            self.injection_potential_at_stations[A] = np.array(self.grabValues(self.injection_potential[A]))

    def getIPModelAndResponse(self, sigma_0, Mn):
        """
        returns the IP model + its responses for given secondary conductivity sigma_0 and chargeaiBlity Mn
        they should be given on integration nodes.
        """
        sigma_oo = sigma_0 + Mn
        txt1 = str(sigma_oo)
        txt2 = str(sigma_0)
        txt3 = str(Mn)
        if getMPIRankWorld() == 0:
            self.logger.debug("sigma_oo = %s" % (txt1))
            self.logger.debug("sigma_0 = %s" % (txt2))
            self.logger.debug("Mn = %s" % (txt3))

        Potentials_0 = {}  # Ws
        Potentials_0_at_Stations = {}
        self.forward_pde.setValue(A=sigma_0 * kronecker(self.forward_pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for iA in self.injection_potential:
            A=self.data.getKeyOfStationNumber(iA)
            self.forward_pde.setValue(X=(self.sigma_background - sigma_0) * grad(self.injection_potential[A]))
            Potentials_0[iA] = self.forward_pde.getSolution()
            Potentials_0_at_Stations[iA] = np.array(self.grabValues(Potentials_0[iA]))
        if getMPIRankWorld() == 0:
            self.logger.debug("%s DC potentials calculated." % len(Potentials_0))
        Potentials_2_at_Stations = {}  # Vs
        Potentials_2 = {}

        self.forward_pde.setValue(A=sigma_oo * kronecker(self.forward_pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for iA in self.injection_potential:
            A = self.data.getKeyOfStationNumber(iA)
            self.forward_pde.setValue(X=Mn * grad(Potentials_0[iA] + self.injection_potential[A]))
            Potentials_2[iA] = self.forward_pde.getSolution()
            Potentials_2_at_Stations[iA] = np.array(self.grabValues(Potentials_2[iA]))

        if getMPIRankWorld() == 0:
            self.logger.debug("%s secondary potentials calculated." % len(Potentials_2))

        return Potentials_0, Potentials_2, Potentials_0_at_Stations, Potentials_2_at_Stations

    def getMisfit(self, sigma_0, Mn, *args):
        """
        return the misfit in potential and chargeaiBlity weighted by weighting_misfit_0 factor
        """
        Potentials_0 = args[0]
        Potentials_2 = args[1]
        Potentials_0_at_Stations = args[2]
        Potentials_2_at_Stations = args[3]

        misfit_2, misfit_0 = 0., 0.
        for iA, iB in self.misfit_2:
            misfit_2 += self.misfit_2[(iA, iB)].getValue(
                Potentials_2_at_Stations[iA] - Potentials_2_at_Stations[iB])  # test against Vs
            misfit_0 += self.misfit_0[(iA, iB)].getValue(
                Potentials_0_at_Stations[iA] - Potentials_0_at_Stations[iB])  # test against Ws
        # print(misfit_2, misfit_0)
        return misfit_2 * self.weightingMisfit_2, misfit_0 * self.weightingMisfit_0

    def getDMisfit(self, sigma_0, Mn, *args):
        """
        returns the derivative of the misfit function with respect to sigma_0 with respect to Mn 
        """
        if len(args) == 0:
            args = self.getIPModelAndResponse(sigma_0, Mn)
        Potentials_0 = args[0]  # Zs
        Potentials_2 = args[1]  # Vs
        Potentials_0_at_Stations = args[2]
        Potentials_2_at_Stations = args[3]

        DMisfitDsigma_0 = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))
        DMisfitDMn = Scalar(0., self.forward_pde.getFunctionSpaceForCoefficient('Y'))

        sigma_oo = (sigma_0 + Mn)
        self.forward_pde.setValue(A=sigma_0 * kronecker(self.forward_pde.getDim()), X=Data(), Y=Data(), y_dirac=Data())
        for iA, iB in self.misfit_2:

            dmis_2 = self.misfit_2[(iA, iB)].getWeightedDifference(
                Potentials_2_at_Stations[iA] - Potentials_2_at_Stations[iB])
            dmis_0 = self.misfit_0[(iA, iB)].getWeightedDifference(
                Potentials_0_at_Stations[iA] - Potentials_0_at_Stations[iB])

            Vs_star_2 = 0
            Q = 0
            for i in range(len(self.misfit_0[(iA, iB)])):
                iM = self.misfit_0[(iA, B)].iMs[i]
                iN = self.misfit_0[(iA, B)].iNs[i]

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
                w1=1.,
                mask_fixed_property = None,
                weighting_misfit_0 = 0.5, pde_tol=1.e-8, stationsFMT="e%s",  logclip=5, **kargs):
        """       
        :domain: pde domain
        :data: survey data, is `fingal.SurveyData`
        :mask_fixed_property: mask of region where Mn and sigma_0 are fixed (set to reference values).
            If not set the bottom of the domain is used.
        :sigma_background: background conductivity, used to calculate the injection potentials
        :sigma_0_ref: reference DC conductivity
        :Mn_ref: reference normalized chargeability (must be positive)
        :param w1: weighting H1 regularization  int grad(m)^2
        :weighting_misfit_0: weighting factor for chargeability part cost function, use 0.5
        :adjVs_ootStationLocationsToElementCenter: moves the station locations to match element centers.
        :stationsFMT: format Vs_ooed to map station keys k to mesh tags stationsFMT%k or None
        :logclip: cliping for p to avoid overflow in conductivity calculation

        """
        super().__init__(domain=domain, data=data, weighting_misfit_0=weighting_misfit_0, sigma_background=sigma_background,
                 pde_tol=pde_tol, stationsFMT=stationsFMT, **kargs)
        self.logclip = logclip
        if mask_fixed_property is None:
            z = self.domain.getX()[2]
            self.mask_fixed_property = whereZero(z - sup(z))
        else:
            self.mask_fixed_property = mask_fixed_property
        self.sigma_0_ref = sigma_0_ref
        self.Mn_ref = Mn_ref



        # regularization
        self.w1 = w1
        self.Hpde = setupERTPDE(self.domain)
        self.Hpde.setValue(A=self.w2 * kronecker(3), q=self.mask_fixed_property)
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

        H1 = self.w2 / 2 * integrate(length(grad(m)) ** 2)
        if getMPIRankWorld() == 0:
            self.logger.debug(
                f'reg: regularization, misfits: ERT, IP =  {H1:e}, {misfit_0:e}, {misfit_2:e}')
        V = H1 + misfit_2 + misfit_0
        return V

    def getGradient(self, m, im, isigma_0, iMn, args2):
        """
        returns the gradient of the cost function. Overwrites `getGradient` of `MeteredCostFunction`
        """
        X = self.w2 * grad(m, where=im.getFunctionSpace())
        DMisfitDsigma_0, DMisfitDMn = self.getDMisfit(isigma_0, iMn, *args2)

        Dsigma_0Dm0 = self.getDsigma0Dm(isigma_0, im[0])
        DMnDm1 = self.getDMnDm(iMn, im[1])
        Y=Data(0., (2,), where=im.getFunctionSpace())
        Y[1] = DMisfitDMn * DMnDm1
        Y[0] = DMisfitDsigma_0 * Dsigma_0Dm0
        return ArithmeticTuple(Y, X)

    def getInverseHessianApproximation(self, r, m, im, isigma_0, iMn, args2, initializeHessian=False):
        """
        returns an approximation of inverse of the Hessian. Overwrites `getInverseHessianApproximation` of `MeteredCostFunction`
        """
        dm=Data(0., (2,), Solution(self.domain))
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
        return integrate(inner(r[0] * m) + inner(r[1], grad(m)))

    def getNorm(self, m):
        """
        returns the norm of property function `m`. Overwrites `getNorm` of `MeteredCostFunction`
        """
        return Lsup(m)