"""
Fingal - some tools 

by l.gross@uq.edu.au, Dec 2020.
"""

from esys.escript import Scalar
from esys.escript import getMPIRankWorld, integrate, hasFeature, Function, kronecker, Data, DiracDeltaFunctions
from esys.escript import inner, length, Lsup, FunctionOnBoundary, inf, sup, whereZero, wherePositive, grad, Solution
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions, LinearPDE
import numpy as np



def makeMaskForOuterSurface(domain, facemask=None, taglist=None):
    """
    returns a mask of elements on the faces where radiation condition is applied.
    if a mask is already given this is returned. If the taglist is given it is used to create this
    otherwise the left, right, back, front and bottom faces are marked.
    """
    if facemask is not None:
        m = facemask
    elif taglist:
        m=Scalar(0., FunctionOnBoundary(domain))
        for t in taglist:
            m.setTaggedValue(t, 1.)
    else:
        X = FunctionOnBoundary(domain).getX()
        x = domain.getX()
        m = wherePositive( whereZero(X[0] - inf(x[0]))+
                           whereZero(X[0] - sup(x[0]))+
                           whereZero(X[1] - inf(x[1]))+
                           whereZero(X[1] - sup(x[1]))+
                           whereZero(X[2] - inf(x[2])))
    return m


def makeWennerArray(numElectrodes=32, id0=0):
    """
    creates a schedule (A,B,M, N) for a Wenner array of length numElectrodes
    """
    schedule = []
    for a in range(1, (numElectrodes + 2) // 3):
        for k in range(0, numElectrodes - 3 * a):
            schedule.append((k + id0, k + 3 * a + id0, k + 1 * a + id0, k + 2 * a + id0))
    return schedule

def makeSchlumbergerArray(numElectrodes=32, id0=0):
    """
    creates a schedule (A,B,M, N) for a Schlumberger array of length numElectrodes
    """
    schedule = []
    for k in range(0, (numElectrodes - 2) // 2):
        for s in range(0, numElectrodes - 2 * k - 3):
            schedule.append((s + id0, s + 2 * k + 3 + id0, s + k + 1  + id0, s + k + 2 + id0))
    return schedule

def makePointSource(source_key, domain, value=1.,stationsFMT=None):
    """

    """
    s = Scalar(0., DiracDeltaFunctions(domain))
    if stationsFMT is None:
        s.setTaggedValue(source_key, value)
    else:
        s.setTaggedValue(stationsFMT % source_key, value)
    return s

def getSourcePotentials(domain, sigma, survey, maskZeroPotential, stationsFMT=None, logger=None):
    """
    return the electric potential for all injections A in the survey using conductivity sigma.
    :domain: PDE domin
    :sigma: conductivity.
    :maskZeroPotential: mask of locations where potentials are zero.
    :return: dictionary of injections A->source_potential assuming conductivity is sigma
    """
    source_potential = {}
    pde = setupERTPDE(domain)
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), q=maskZeroPotential)
    for A in survey.getStationNumeration() : # getListOfInjectionStations():
        iA=survey.getStationNumber(A)
        pde.setValue(y_dirac=makePointSource(A, pde.getDomain(),stationsFMT=stationsFMT))
        source_potential[iA] = pde.getSolution()
        #if logger:
        #      logger.debug(f"processing source {iA} key = {A}: {source_potential[iA]}")
        assert Lsup(source_potential[iA]) > 0, "Zero potential for injection %s" % A
    if logger:
       logger.debug(f"{len(source_potential)} source potentials calculated.")
    return source_potential



def getSourcePotentialsWithRobinCondition(domain, sigma, survey, sigma_faces=None, maskZeroPotential=None, stationsFMT=None, logger=None):
    """
    return the electric potential for all injections A in the survey using conductivity sigma.
    :sigma: conductivity. Needs to a constant if sigma_source is not present.
    :sigma_faces: conductivity at faces. If None sigma is used.
    :maskOuterFaces:
    :return: dictionary of injections A->injection_potentials assuming internal and surface conductivity is sigma and
             sigma_surface
    """
    if not sigma_faces:
        sigma_faces = sigma
    mm = makeMaskForOuterSurface(domain, facemask=maskOuterFaces)
    n = domain.getNormal() * mm

    x = n.getX()
    source_potential = {}
    pde = setupERTPDE(domain)
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
    for A in survey.getListOfInjectionStations():
        iA=survey.getStationNumber(A)

        pde.setValue(y_dirac=makePointSource(A, pde.getDomain(),stationsFMT=stationsFMT))
        #pde.setValue(y=1)
        # ---
        xA = survey.getStationLocationByKey(A)
        r = x - xA
        pde.setValue(d=sigma_faces * inner(r, n) / length(r) ** 2)  # doi:10.1190/1.1440975
        source_potential[iA] = pde.getSolution()
        #if logger:
        #      logger.debug(f"processing source {iA} key = {A}: {source_potential[iA]}")
        assert Lsup(source_potential[iA]) > 0, "Zero potential for injection %s" % A
    if logger:
       logger.debug(f"{len(source_potential)} source potentials calculated.")
    return source_potential

def getAdditivePotentials(pde, sigma, schedule, sigma_stations=None, source_potential={},
                          sigma_src=1., sigma_src_stations=None, logger=None):
    """
    calculates the extra/additive potentials V_A for given sigma and source potentials for sigma_src

    :param pde: PDE used to get the secondary potential. Coefficients are altered.
    :param sigma: electrical conductivity
    :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
    :type schedule: `SurveyData`
    :param sigma_stations:electric conductivity at_stations
    :type sigma_stations:  `dict` of `float`
    :param source_potential: potentionals for a sources with assumed conductivity sigma_src
    :type source_potential: `dict` of `Scalar`
    :param sigma_src: conductivity used to calculate source potentials
    :param sigma_src_stations: conductivity for source potentials  at_stations. If not set  'sigma_src' is used.
    :param maskZeroPotential: mask where additional potential is to zero.
    :type maskZeroPotential: `Scalar` with `Solution` attribute
    :return:  secondary/incremental potential to source potential due to conductivity sigma and the associated rescaling
                factor to adjust conductivity at source.
    """
    if sigma_src_stations is None:
        sigma_src_stations = {iA: sigma_src for iA in source_potential}

    potential = {}
    src_potential_scale = {}
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data())
    for iA in source_potential:
        if sigma_stations is None:
            alpha_A = 1.
        elif isinstance(sigma_src_stations, dict) or isinstance(sigma_src_stations, list):
            alpha_A = sigma_src_stations[iA] / sigma_stations[iA]
        else:
            alpha_A = sigma_src_stations / sigma_stations[iA]
        pde.setValue(X=(sigma_src - sigma * alpha_A) * grad(source_potential[iA]))
        potential[iA] = pde.getSolution()
        src_potential_scale[iA] = alpha_A
    if logger:
        logger.info(f"{len(potential)} additive DC potentials calculated.")
    return potential, src_potential_scale


def getAdditivePotentialsWithRobinCondition(pde, sigma, sigma_faces, schedule, sigma_stations=None, source_potential={},
                          sigma_src=1., sigma_src_faces=None, sigma_src_stations=None, mask_faces=None,
                          logger=None):
    """
    calculates the extra/additive potentials V_A for given sigma and source potentials for sigma_src

    :param pde: PDE used to get the secondary potential. Coefficients are altered.
    :param sigma: electrical conductivity
    :param sigma_faces: electrical conductivity at faces
    :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
    :type schedule: `SurveyData`
    :param sigma_stations:electric conductivity at_stations
    :type sigma_stations:  `dict` of `float`
    :param source_potential: potentionals for a sources with assumed conductivity sigma_src
    :type source_potential: `dict` of `Scalar`
    :param sigma_src: conductivity used to calculate source potentials
    :param sigma_src_faces: conductivity for source potentials  at surfaces. If not set  'sigma_src' is used.
    :param sigma_src_stations: conductivity for source potentials  at_stations. If not set  'sigma_src' is used.
    :param mask_faces: mask of surface elements to apply `radiation` conditions.
        If not set to bottom, front, back, left, right elements.
    :type mask_faces: `None` or `Scalar` with `FunctionOnBoundary` attribute
    :return:  secondar/incremental potential to source potential due to conductivity sigma
    """
    if mask_faces is None:
        mask_faces = makeMaskForOuterSurface(pde.getDomain())
    if sigma_src_faces is None:
        sigma_src_faces = sigma_src
    if sigma_src_stations is None:
        sigma_src_stations = {iA: sigma_src for iA in source_potential}
    n = pde.getDomain().getNormal() * mask_faces
    x_bc = FunctionOnBoundary(pde.getDomain()).getX()


    potential = {}
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), d=Data(), X=Data(), Y=Data(), y=Data())
    for iA in source_potential:
        if sigma_stations is None:
            alpha_A = 1.
        elif isinstance(sigma_src_stations, dict):
            alpha_A = sigma_src_stations[iA] / sigma_stations[iA]
        else:
            alpha_A = sigma_src_stations / sigma_stations[iA]
        xA = schedule.getStationLocationByNumber(iA)
        r = x_bc - xA
        fA = inner(r, n) / length(r) ** 2
        pde.setValue(d=sigma_faces * fA, y=(sigma_src_faces - sigma_faces * alpha_A) * fA * source_potential[iA])
        pde.setValue(X=(sigma_src - sigma * alpha_A) * grad(source_potential[iA]))
        potential[iA] = pde.getSolution() + (alpha_A - 1.) * source_potential[iA]
    if logger:
        logger.info(f"{len(potential)} additive DC potentials calculated.")
    return potential

def getSecondaryPotentials(pde, sigma_oo, M_n, source_potential={}, src_potential_scale_DC = {}, additive_potential_DC ={}, logger=None):
    """
    calculates the extra/additive potentials V_A for given sigma and source potentials for sigma_src

    :param pde: PDE used to get the secondary potential. Coefficients are altered.
    :param sigma_oo: electrical conductivity (for high frequencies)
    :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
    :type schedule: `SurveyData`
    :param M_n: normalized chargeability
    :type M_n:  `dict` of `Scalar`
    :param source_potential: potentials for a sources with assumed conductivity
    :type source_potential: `dict` of `Scalar`
    :param additive_potential_DC: additive potentials for a sources for sigma_0 = sigma_oo - M_n
    :type additive_potential_DC: `dict` of `Scalar`
    :type mask_faces: `None` or `Scalar` with `FunctionOnBoundary` attribute
    :return:  secondary/over-voltage potential
    """
    potential = {}
    pde.setValue(A=sigma_oo * kronecker(3), y_dirac=Data(), X=Data(), Y=Data())
    for iA in source_potential:
        VA = src_potential_scale_DC[iA]  * source_potential[iA] + additive_potential_DC[iA]
        pde.setValue(X=M_n * grad(VA))
        potential[iA] = pde.getSolution()
    logger.info(f"{len(potential)} secondary IP potentials calculated.")
    return potential

def makeZZArray(numElectrodes=32, id0=0):
    """
    creates a schedule (A,B,M, N) for a Data0 array of length numElectrodes (full monty)
    """
    schedule = []
    for a in range(numElectrodes):
        for b in range(a + 1, numElectrodes):
            for m in range(numElectrodes):
                for n in range(m + 1, numElectrodes):
                    if set([a, b]).isdisjoint([m, n]):
                        schedule.append((a + id0, b + id0, m + id0, n + id0))
    return schedule


def FindNearestElectrode(x, y, z, electrodes={}):
    """
    finds the nearest electrode in the dictionary electrodes
    (x,y,z) is the coordinates of a proposed electrode
    :param x: x coordinate of electrode
    :param y: y coordinate of electrode
    :param z: x coordinate of electrode
    :param electrodes: dictionary of defined electrodes: electrodes[ide]=X=(xe, ye, ze)
    :return: id of nearest electrode and distance and distance to 
    
    """
    distmin = 1e88
    idemin = None
    for ide in electrodes:
        X = electrodes[ide]
        dist = ((X[0] - x) ** 2 + (X[1] - y) ** 2 + (X[2] - z) ** 2) ** 0.5
        if dist < distmin:
            distmin = dist
            idemin = ide
    if not idemin is None:
        return int(idemin), distmin
    else:
        return idemin, distmin

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
        data_log = np.log(abs(data)+self.EPS)
        super().__init__(iMs, data_log, iNs, weightings, injections, **kwargs)
    def getValue(self, u):
        res = self.data-np.log(abs(self.getDifference(u))+self.EPS)
        dd = abs(res) ** 2 * self.weightings
        return 0.5 * sum(dd)

    def getDerivative(self, u):
        nn=self.getDifference(u)
        res=self.data-np.log(abs(nn)+self.EPS)
        return -res * self.weightings/nn

def setupPDESystem(domain, numEquations=1, symmetric=True, tolerance=1e-8):
    pde = LinearPDE(domain, numEquations=numEquations, numSolutions=numEquations, isComplex=False)
    if symmetric:
        pde.setSymmetryOn()
    optionsG = pde.getSolverOptions()
    optionsG.setSolverMethod(SolverOptions.PCG)
    optionsG.setTolerance(tolerance)

    if hasFeature('trilinos'):
        optionsG.setPackage(SolverOptions.TRILINOS)
        optionsG.setPreconditioner(SolverOptions.AMG)
        optionsG.setTrilinosParameter("verbosity", "none")
        optionsG.setTrilinosParameter("number of equations", numEquations)
        optionsG.setTrilinosParameter("problem: symmetric", symmetric)

    return pde
def setupERTPDE(domain, tolerance=1e-8, poisson=True, debug=0):
    """
    used to setup all PDEs fro ERT related inversion. If available TRILINOS is usered.
    
    :param domain: domain of the PDE
    :type domain: `esys.escript.AbstractDomain`
    :param tolerance: solver tolerance
    :type tolerance: `float` 
    :param poisson: if True TRILINOS settings fro Poisson problems is set.
    :type poisson: `bool`
    :return: linear, scalar PDE with real coefficients.
    :rtype: `esys.escript.linearPDEs.LinearPDE`
    

    """
    pde = LinearSinglePDE(domain, isComplex=False)
    pde.setSymmetryOn()
    optionsG = pde.getSolverOptions()
    optionsG.setSolverMethod(SolverOptions.PCG)
    optionsG.setTolerance(tolerance)
    #optionsG.setSolverMethod(SolverOptions.DIRECT)
    if hasFeature('trilinos'):
        optionsG.setPackage(SolverOptions.TRILINOS)
        optionsG.setPreconditioner(SolverOptions.AMG)
        if poisson:
            optionsG.setTrilinosParameter("problem:type", "Poisson-3D")
        optionsG.setTrilinosParameter("verbosity", "none")
        optionsG.setTrilinosParameter("number of equations", 1)
        optionsG.setTrilinosParameter("problem: symmetric", True)
    # print("TO DO")
    #optionsG.setSolverMethod(SolverOptions.DIRECT)

    return pde


def getR2(y, y_fitted, chi=None):
    """
        calculates the coefficient of determination R^2 for  `y_fitted` as prediction for `y` over a region marked by chi>0 defined by
        
        R^2=1 - S_res/S_tot
        
        with S_res=int(chi*(y-y_fitted*1)**2,  S_tot=int(chi*(y-m(y)*1)**2), m(y)=int(chi*y)/int(chi)
        
        If R^2=1 then `y_fitted` is predicts `y` exactly. If R^2 then `y_fitted` does not make a better prediction than the mean. 
        
        :param y: target distribution
        :type y: `esys.escript.Scalar`
        :param y_fitted: fitted distribution
        :type y_fitted: `esys.escript.Scalar`
        :param chi: marker/weighting for region of interest
        :type chi: `esys.escript.Scalar` or None
        :rtype: `float`
        """
    if chi is None:
        chi = Scalar(1., Function(y_fitted.getFunctionSpace().getDomain()))

    ybar = integrate(chi * y) / integrate(chi)
    S_res = integrate(chi * (y - y_fitted) ** 2)
    S_tot = integrate(chi * (y - ybar) ** 2)
    if S_tot > 0:
        R2 = 1 - S_res / S_tot
    else:
        if S_res > 0:
            R2 = 0.
        else:
            R2 = 1.

    return R2

class InterpolatorWithExtension(object):
    """
    An interpolator of 2D cloud data. The data are first interpolated to a rectangular grid
    using nearest naighbour interpolation. If markExtrapolation is values obtained in the interpolation
    from the rectangular grid using extrapolated values are marked. The grid spacing is determined by the
    smallest distance of the location of data point.

    :param locations: 2D array of locations where data aer given
    :param data: data values at locations
    :param xmin: minimal X coordinate of interpolation domain
    :param xmax: maximal X coordinate of interpolation domain
    :param ymin: minimal Y coordinate of interpolation domain
    :param ymax: maximal Y coordinate of interpolation domain
    :param markExtrapolation: If values obtained by extrapolation are marked.
    """
    def __init__(self, locations, data, xmin, xmax, ymin, ymax, markExtrapolation=False):
        from scipy.interpolate import griddata, NearestNDInterpolator, RegularGridInterpolator
        import numpy as np
        self.markExtrapolation = markExtrapolation
        X, Y = locations[0], locations[1]
        dx = 1e99
        TOL = abs(X).max() * 1e-8
        for i in range(len(X) - 1):
            p = abs(X[i] - X[i + 1:]) + abs(Y[i] - Y[i + 1:])
            dx = min(p[p > 0.].min(), dx)

        nx = int((xmax - xmin) / dx + 0.5)
        ny = int((ymax - ymin) / dx + 0.5)
        self.grid = (nx, ny)
        self.dx = dx

        print(f"Interpolation grid is {nx} x {ny} of spacing {dx}.")
        X = np.linspace(xmin, xmax, num=nx)
        Y = np.linspace(ymin, ymax, num=ny)
        grid_x, grid_y = np.meshgrid(X, Y, indexing='ij')
        if self.markExtrapolation:
            z_grid = griddata(locations, data, (grid_x, grid_y), method='nearest')
            test_grid = griddata(locations, data, (grid_x, grid_y), method='linear')
            self.interpolator_to_test_extrapolation = RegularGridInterpolator((X, Y), test_grid)
        else:
            z_grid = griddata(locations, data, (grid_x, grid_y), method='nearest')
        self.interpolator = RegularGridInterpolator((X, Y), z_grid)

    def __call__(self, x, y):
        """
        interpolate to location (x,y) within the interpolation range.
        the value together with index 0 or 1 for 'extrapolated' or 'genuine' is returned.
        """
        if self.markExtrapolation:
            v = self.interpolator_to_test_extrapolation((x, y))
            i = 1
            if np.isnan(v):
                i = 0
                v = self.interpolator((x, y))
            return v, i
        else:
            return self.interpolator((x, y)), 1


def mapToDomain2D(target, interpolator, where=None):
    """
    this the interpolator(x,y) to fill ``Data`` object `target`.
    `target` annd `target_interpolated` are returned where
     values in `target` with `target_interpolated` >0 obtained by interpolation.
    If `where` is given only location with `where` value > 0 are set.
    """
    X = target.getX()
    target.expand()
    target_interpolated = Scalar(0., target.getFunctionSpace())
    target_interpolated.expand()
    if where is None:
        for p in range(X.getNumberOfDataPoints()):
            x = X.getTupleForDataPoint(p)
            v, ind = interpolator(x[0], x[1])
            target.setValueOfDataPoint(p, v)
            target_interpolated.setValueOfDataPoint(p, ind)
    else:
        assert target.getFunctionSpace() == where.getFunctionSpace()
        for p in range(X.getNumberOfDataPoints()):
            if where.getTupleForDataPoint(p)[0] > 0:
                x = X.getTupleForDataPoint(p)
                v, ind = interpolator(x[0], x[1])
                target.setValueOfDataPoint(p, v)
                target_interpolated.setValueOfDataPoint(p, ind)

    return target, target_interpolated

def createBackgroundTemperature(domain, T=20., mask_fixed_temperature=None, gradT=0.03, mask_top_surface=None):
    """
    create a background temperature from a temperature `T` at locations marked by `mask_fixed_temperature`.
    On the top and bottom a temperature gradient `gradT` is assumed which might be overwritten by the set temperature.

    :param domain: domain
    :param T: given surface temperature
    :param mask_fixed_temperature: mask of locations where temperature is set (on `Solution`).
                                 If not set the (flat) top surface is used..
    :param gradT: temperature gradient
    :param mask_top_surface: mask top surface (on FunctionOnBoundary). If not set the (flat) top surface is used.
    :return: a temperature profile.
    """
    z=domain.getX()[2]
    z_inf=inf(z)
    z_sup=sup(z)
    Z=FunctionOnBoundary(domain).getX()[2]
    if  mask_top_surface is None:
        mask_top_surface=whereZero(Z-z_sup)
    mask_bottom_surface=whereZero(Z-z_inf)
    Q=-mask_top_surface*gradT+mask_bottom_surface*gradT
    if  mask_fixed_temperature is None:
        mask_fixed_temperature=whereZero(z-z_sup)
    pde=LinearSinglePDE(domain, isComplex=False)
    pde.setSymmetryOn()
    optionsG=pde.getSolverOptions()
    #optionsG.setTolerance(1e-10)
    pde.setValue(A=kronecker(3), y=Q, r=T, q=mask_fixed_temperature)
    Tnew=pde.getSolution()
    return Tnew

from esys.escript.pdetools import Locator
from esys.escript import interpolate

class PoissonEquationZeroMean(object):
    """
    solves the Poisson equation for with a secondary normalization condition (h linear):
            -u_,ii = -div(X) + Y
            n_i*u_,i = y
            h(u) = 0

    In weak form (with Lagrange Multiplier) Find u and t with: for all v:

        a(v, u) + t * h(v) = f(v)
        h(u)=hu

    with: a(v, u) = Int( A_ij v_,i u_,i)
        f(v) = Int( X_i * v_,i  + Y * v ) + int(y * v)
        h(v) = Int( Xh_i * v_,i  + Yh * v ) + int(yh * v)

    Int = Volume integral, int = surface integral

    Solution method:
        pick point X_0.
        assume for a FEM basis function phi_0 we can write: v = v0 + v(X0) * phi0 with  v0(X0)=0
        We have a(v+c, u) = a(u,, v+c ) =0 for any constant c

        solve uf(X0)=0, for all v0 with  v0(X0)=0: a(v0, uf) =  f(v0)
        solve uh(X0)=0, for all v0 with  v0(X0)=0: a(v0, uh) =  h(v0)

        then set:

                u = uf - t * uh + c

        Then by construction of u0 and uh: for all v0 with v0(X_0)=0
                a(v0, u) + t * h(v0)  =  a(v0, uf - t * uh) + t * h(v0) = f(v0)
        also
                a(phi0, u) + t * h(phi0)  =  a(phi0, uf) - t * a(phi_0, uh) + t * h(phi0) = f(phi0)
        if we set:
            t = (f(phi0)-a(phi0, uf))/(h(phi0)-a(phi0, uh)) = f1/h1

        we have h(u) = hu if we set

            c = (hu -h(uf - t * uh) ) /h1

        Notes:  h1 = h(1) = Int(Yh) + int(yh)
                h(phi0)-a(phi0, uh) = h(1) = h1 (as phi0 = 1 - v0 with v0(X_0)=0)
                f(phi0)-a(phi0, uf) = f(1) = f1 (again as phi0 = 1 - v0 with v0(X_0)=0)

    """
    def __init__(self, domain, X0=None, A = None, rtol = 1e-8, isComplex = False):
        self.domain=domain
        numDim = domain.getDim()
        self.pde = LinearSinglePDE(domain, isComplex=isComplex)
        self.pde.setSymmetryOn()
        optionsG = self.pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.PCG)
        optionsG.setTolerance(rtol)
        # optionsG.setSolverMethod(SolverOptions.DIRECT)
        if hasFeature('trilinos'):
            optionsG.setPackage(SolverOptions.TRILINOS)
            optionsG.setPreconditioner(SolverOptions.AMG)
            if numDim == 3:
                optionsG.setTrilinosParameter("problem:type", "Poisson-3D")
            else:
                optionsG.setTrilinosParameter("problem:type", "Poisson-2D")
            optionsG.setTrilinosParameter("verbosity", "none")
            optionsG.setTrilinosParameter("number of equations", 1)
            optionsG.setTrilinosParameter("problem: symmetric", True)
        x = Solution(domain).getX()
        if not X0:
            X0 = tuple( [ inf(x[i]) for i in range(numDim) ] )
        self.X0  = Locator(Solution(domain),  X0).getX()
        q=whereZero(length(x-self.X0))
        if A is None :
            self.pde.setValue(A = kronecker(numDim), q=q)
        else:
            self.pde.setValue(A=A, q=q)

    def getX0(self):
        """
        return the used location of X0
        """
        return self.__loc.getX()
    def getPDE(self):
        """
        return the PDE to be solved
        """
        return self.pde

    def setSecondaryCondition(self, Yh=None, Xh=None, yh=None ):

        self.h1 = 0
        if Yh is None:
            self.Yh = Scalar(1, Function(self.domain))
        else:
            self.Yh = interpolate(Yh, Function(self.domain))
        self.h1+= integrate(self.Yh)
        if Xh is None:
            self.Xh = Data()
        else:
            self.Xh = interpolate(Xh, Function(self.domain))
        if yh is None:
            self.yh = Data()
        else:
            self.yh = interpolate(yh, FunctionOnBoundary(self.domain))
            self.h1+= integrate(self.yh)
        assert abs(self.h1) > 0 , "h1=h(1) must be non-zero."
        # ... recover uh
        self.pde.setValue(X=self.Xh, Y=self.Yh, y=self.yh)
        self.uh=self.pde.getSolution()
        return self

    def getSolution(self, X=Data(), Y=Data(), y=Data(), hu=0):
        # ... recover uh
        self.pde.setValue(X=X, Y=Y, y=y)
        uf = self.pde.getSolution()
        f1=0
        if not Y.isEmpty():
            f1+=integrate(self.pde.getCoefficient("Y"))
        if not y.isEmpty():
            f1 += integrate(self.pde.getCoefficient("y"))
        t = f1/self.h1
        u = uf - t * self.uh
        hu1=0
        if not self.Xh.isEmpty():
            hu1+=integrate(inner(self.Xh, grad(u, self.Xh.getFunctionSpace() )))
        if not self.Yh.isEmpty():
            hu1+=integrate(self.Yh * u )
        if not self.yh.isEmpty():
            hu1+=integrate(self.yh * u )
        c = (hu - hu1)/self.h1
        u+=c
        return u