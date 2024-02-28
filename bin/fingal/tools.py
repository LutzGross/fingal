"""
Fingal - some tools 

by l.gross@uq.edu.au, Dec 2020.
"""

from esys.escript import Scalar, getMPIRankWorld, integrate, hasFeature, Function, kronecker, Data, DiracDeltaFunctions, \
    inner, length, Lsup, FunctionOnBoundary, inf, sup, whereZero, wherePositive, grad
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions


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


def makePointSource(source_key, domain, value=1., stationsFMT=None):
    """

    """
    s = Scalar(0., DiracDeltaFunctions(domain))
    if stationsFMT is None:
        s.setTaggedValue(source_key, value)
    else:
        s.setTaggedValue(stationsFMT % source_key, value)
    return s


def getSourcePotentials(domain, sigma, survey, sigma_surface=None, mask_outer_faces=None, stationsFMT=None):
    """
    return the electric potential for all injections A in the survey using conductivity sigma.
    :sigma: conductivity. Needs to a constant if sigma_source is not present.
    :sigma_surface: surface condcutivity. If None sigma is used.
    :mask_outer_faces
    :return: dictionary of injections A->injection_potentials assuming internal and surface conductivity is sigma and
             sigma_surface
    """
    if not sigma_surface:
        sigma_surface = sigma
    mm = makeMaskForOuterSurface(domain, facemask=mask_outer_faces)
    n = domain.getNormal() * mm
    x = n.getX()

    source_potential = {}
    pde = setupERTPDE(domain)
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
    for A in survey.getListOfInjectionStations():
        iA=survey.getStationNumber(A)
        pde.setValue(y_dirac=makePointSource(A, pde.getDomain(), stationsFMT=stationsFMT))
        # ---
        xA = survey.getStationLocationByKey(A)
        r = x - xA
        pde.setValue(d=sigma * inner(r, n) / length(r) ** 2)  # doi:10.1190/1.1440975
        source_potential[iA] = pde.getSolution()
        print("processing source ", iA, " key = ", A, source_potential[iA])
        assert Lsup(source_potential[iA]) > 0, "Zero potential for injection %s" % A
    return source_potential


def getSecondaryPotentials(pde, sigma, sigma_at_face, schedule, sigma_at_station=None, source_potential={},
                           sigma_src=1., sigma_src_at_face=None, sigma_src_at_station=None, mask_faces=None):
    """
    calculates the extra/secondary potentials V_A for given sigma and source potentials for sigma_src

    :param pde: PDE used to get the secondary potential. Coefficients are altered.
    :param sigma: electrical conductivity
    :param sigma_at_face: electrical conductivity at faces
    :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
    :type schedule: `SurveyData`
    :param sigma_at_station:electric conductivity at stations
    :type sigma_at_station:  `dict` of `float`
    :param source_potential: potentionals for a sources with assumed conductivity sigma_src
    :type source_potential: `dict` of `Scalar`
    :param sigma_src: conductivity used to calculate source potentials
    :param sigma_src_at_face: conductivity for source potentials  at surfaces. If not set  'sigma_src' is used.
    :param sigma_src_at_stations: conductivity for source potentials  at stations. If not set  'sigma_src' is used.
    :param mask_faces: mask of surface elements to apply `radiation` conditions.
        If not set to bottom, front, back, left, right elements.
    :type mask_faces: `None` or `Scalar` with `FunctionOnBoundary` attribute
    :return:  secondar/incremental potential to source potential due to conductivity sigma
    """
    if mask_faces is None:
        mask_faces = makeMaskForOuterSurface(pde.getDomain())
    if sigma_src_at_face is None:
        sigma_src_at_face = sigma_src
    if sigma_src_at_station is None:
        sigma_src_at_station = {iA: sigma_src for iA in source_potential}
    n = pde.getDomain().getNormal() * mask_faces
    x_bc = FunctionOnBoundary(pde.getDomain()).getX()

    potential = {}
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
    for iA in source_potential:
        if sigma_at_station is None:
            alpha_A = 1.
        elif isinstance(sigma_src_at_station, dict):
            alpha_A = sigma_src_at_station[iA] / sigma_at_station[iA]
        else:
            alpha_A = sigma_src_at_station / sigma_at_station[iA]
        xA = schedule.getStationLocationByNumber(iA)
        r = x_bc - xA
        fA = inner(r, n) / length(r) ** 2
        pde.setValue(d=sigma_at_face * fA, y=(sigma_src_at_face - sigma_at_face * alpha_A) * fA * source_potential[iA])
        pde.setValue(X=(sigma_src - sigma * alpha_A) * grad(source_potential[iA]))
        potential[iA] = pde.getSolution() + (alpha_A - 1.) * source_potential[iA]
        print("processing station number ", iA, ": alpha=", alpha_A, "sol. V ", pde.getSolution(), " DV =", potential[iA])

    return potential

def makeZZArray(numElectrodes=32, id0=0):
    """
    creates a schedule (A,B,M, N) for a ZZ array of length numElectrodes (full monty)
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
    if hasFeature('trilinos'):
        if debug and getMPIRankWorld() == 0:
            print("TRILINOS solver used.")
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
