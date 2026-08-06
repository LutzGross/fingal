"""
Smoothed continuous damage model for 3-D brittle materials.

Implements the local continuum damage constitutive model of

    S. Mondal, L.M. Olsen-Kettle, L. Gross,
    "Regularization of continuum damage mechanics models for 3-D brittle
     materials using implicit gradient enhancement",
    Computers and Geotechnics 122 (2020) 103505.
    https://doi.org/10.1016/j.compgeo.2020.103505

The model is expressed in terms of `esys.escript` `Data` objects so that it
works unchanged from a single finite element up to large unstructured meshes.
The non-local (implicit gradient) extension replaces the local equivalent
strain fed into `update()` by the solution of the Helmholtz equation
`ebar - c*laplace(ebar) = etilde`; on a single element (Neumann boundary,
uniform field) local and non-local strain coincide, so this class already
captures the full single-element behaviour.

by fingal, 2026.
"""

import logging
import numpy as np
from esys.escript import (Function, Solution, Scalar, Vector, Data, kronecker,
                          trace, inner, sqrt, clip, maximum, minimum, symmetric,
                          grad, Lsup, sup, inf, integrate, interpolate)
from esys.escript.linearPDEs import LinearSinglePDE, LameEquation, SolverOptions

__all__ = ['SmoothDamageModel']


def _asField(value, fs):
    """
    returns `value` as a `Data` object on the function space `fs`, accepting
    either a scalar (uniform material) or a `Data` object (heterogeneous
    material, e.g. a layered specimen).
    """
    if isinstance(value, Data):
        return interpolate(value, fs)
    return Scalar(float(value), fs)


def _inf(value):
    """
    minimum of `value` over the domain; passes scalars through unchanged. Used
    to validate material parameters that may be given as `Data`.
    """
    return inf(value) if isinstance(value, Data) else float(value)


def _sup(value):
    """
    maximum of `value` over the domain; passes scalars through unchanged.
    """
    return sup(value) if isinstance(value, Data) else float(value)


def _rootIncreasing(phi, x_scale, tol, itmax=50):
    """
    finds the root of a monotonically increasing scalar function `phi` near 0 by
    bracketing (expanding a step of magnitude `|x_scale|` in the direction that
    reduces `phi`) followed by bisection to the absolute tolerance `tol`. Returns
    the (possibly negative) root; `phi` is evaluated by field operations, so calls
    are kept modest. Used for the dissipation-arc-length step size.
    """
    f0 = phi(0.)
    if abs(f0) <= tol:
        return 0.
    dx = abs(x_scale) if x_scale != 0. else 1e-6
    if f0 < 0.:                                   # root at positive argument
        lo, flo, hi = 0., f0, dx
        for _ in range(60):
            fhi = phi(hi)
            if fhi > 0.:
                break
            lo, hi = hi, 2. * hi
        else:
            return hi
    else:                                         # root at negative argument
        hi, fhi, lo = 0., f0, -dx
        for _ in range(60):
            flo = phi(lo)
            if flo < 0.:
                break
            hi, lo = lo, 2. * lo
        else:
            return lo
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fm = phi(mid)
        if abs(fm) <= tol:
            return mid
        if fm < 0.:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def isotropicStiffnessTensor(E, nu):
    """
    returns the isotropic linear-elastic 4th order stiffness tensor C_ijkl
    (`numpy` array of shape (3,3,3,3)) for Young's modulus `E` and Poisson
    ratio `nu`, using the Lame parameters

        lam = E*nu/((1+nu)*(1-2*nu)),  mu = E/(2*(1+nu)).
    """
    lam = E * nu / ((1. + nu) * (1. - 2. * nu))
    mu = E / (2. * (1. + nu))
    d = np.eye(3)
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i, j, k, l] = (lam * d[i, j] * d[k, l]
                                     + mu * (d[i, k] * d[j, l] + d[i, l] * d[j, k]))
    return C


class SmoothDamageModel(object):
    """
    Local continuum damage model with modified von Mises equivalent strain and
    modified power-law damage evolution.

    The damage variable D in [0,1] scales the elastic stiffness as
    C = (1-D)*C0. It is an explicit function of the monotonically increasing
    history variable kappa which tracks the largest equivalent strain seen so
    far (never below the damage threshold kappa0).

    Default parameters follow Table 1 (damage law) and Table 2 (elastic
    properties, strength ratio) of Mondal et al. (2020).

    All material parameters (`E0`, `nu`, `kappa0`, `kappa_c`, `alpha`, `beta`,
    `gamma`, `sigma_t`, ...) may be given either as scalars for a uniform
    material or as `Data` objects on `Function(domain)` for a heterogeneous one
    (e.g. the layered roof/coal/floor specimen of
    `examples/SimpleCompression`). Note that heterogeneous parameters must be
    built on the same domain that is later passed to `initialize`.
    """
    def __init__(self, E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                 alpha=1.6, beta=0.01, gamma=841. / 250., sigma_t=None,
                 min_stiffness_ratio=1.e-3, localization_length=None,
                 porosity0=0.0, biot=None, crack_porosity=0.0):
        """
        :param E0: undamaged Young's modulus [Pa]
        :param nu: Poisson ratio
        :param kappa0: threshold (damage-initiation) equivalent strain
        :param kappa_c: ultimate equivalent strain (D=1 for kappa>=kappa_c)
        :param alpha: softening parameter (final softening stage)
        :param beta: damage growth influence parameter (initial growth rate)
        :param gamma: ratio of compressive to tensile strength sigma_c/sigma_t
        :param sigma_t: tensile strength [Pa]; if given, kappa0 is overwritten
                        by kappa0 = sigma_t/E0 (Eq. 11).
        :param min_stiffness_ratio: residual stiffness fraction (1-D floored at
                        this value in `getStiffnessTensor`) to keep the
                        equilibrium system non-singular as D -> 1.
        :param localization_length: localization length scale l [m] of the
                        implicit-gradient (non-local) enhancement. If set, the
                        equivalent strain is smoothed by the Helmholtz equation
                        ebar - c*laplace(ebar) = etilde with c = l**2/2 before
                        it drives the damage. If `None` (default) the local
                        model is used (ebar = etilde).
        :param porosity0: initial porosity phi0 (auxiliary field, see
                        `getPorosity`). Default 0 leaves porosity tracking off.
        :param biot: intact Biot coefficient b0 = 1 - K/Ks (drained skeleton bulk
                        modulus over solid-grain bulk modulus) of the poroelastic
                        volumetric porosity term. Defaults to 1 - porosity0.
        :param crack_porosity: coefficient `a` of the (dilatant, damage-driven)
                        crack porosity a*D added in `getPorosity`.
        """
        if sigma_t is not None:
            kappa0 = sigma_t / E0
        assert 0. <= _inf(nu) and _sup(nu) < 0.5, \
            "Poisson ratio must be in [0, 0.5)."
        assert _inf(kappa0) > 0., "need kappa0 > 0."
        assert _inf(kappa_c - kappa0) > 0., "need kappa_c > kappa0 everywhere."
        assert _inf(gamma) >= 1., \
            "compressive strength must not be below tensile strength."
        self.E0 = E0
        self.nu = nu
        self.kappa0 = kappa0
        self.kappa_c = kappa_c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.min_stiffness_ratio = min_stiffness_ratio
        self.localization_length = localization_length
        self.porosity0 = porosity0
        self.biot = (1. - porosity0) if biot is None else biot   # intact Biot b0
        self.crack_porosity = crack_porosity
        self.porosity = None                     # auxiliary porosity state field
        self.logger = logging.getLogger('fingal.SmoothDamageModel')
        # isotropic Lame parameters of the undamaged material; the full 4th
        # order stiffness tensor is `isotropicStiffnessTensor(E0, nu)`.
        self.lam = E0 * nu / ((1. + nu) * (1. - 2. * nu))
        self.mu = E0 / (2. * (1. + nu))
        self.kappa = None
        self.dim = None                          # set by `initialize`
        self.elasticity = None
        self.helmholtz = None
        # accumulated total displacement and total Cauchy stress (state carried
        # between load steps by the incremental scheme); set in `initialize`.
        self.u = None
        self.stress = None
        # diagnostics of the most recent `solveLoadStep`.
        self.converged = False
        self.last_change = None

    def initialize(self, domain, elasticity_solver=SolverOptions.PCG):
        """
        prepares the model on `domain`: sets the history variable kappa to the
        threshold value kappa0 and the damage D to zero on the element
        (`Function`) function space, builds the `LameEquation` for the
        elasticity problem (available as `self.elasticity`) and, if a
        localization length is set, the Helmholtz PDE for the non-local
        smoothing. Returns self.

        The caller sets the constraints (`q`, `r`) on `self.elasticity` for the
        loading; the damaged Lame parameters are (re)set by `solveLoadStep`. The
        accumulated total displacement `self.u` and total stress `self.stress`
        (the state carried by the incremental solve) are reset to zero.

        :param elasticity_solver: solver method for the elasticity PDE. An
                        iterative method (e.g. PCG, the default) is preconditioned
                        with algebraic multigrid (AMG); use DIRECT for tiny
                        problems (e.g. the single-element test).
        """
        dim = domain.getDim()
        self.dim = dim
        # `_asField` accepts scalar or `Data` parameters (heterogeneous material)
        self.kappa = _asField(self.kappa0, Function(domain))
        self.D = Scalar(0., Function(domain))
        self.porosity = _asField(self.porosity0, Function(domain))  # auxiliary
        self.u = Vector(0., Solution(domain))                # total displacement
        self.stress = Data(0., (dim, dim), Function(domain))  # total Cauchy stress
        # elasticity problem  -(sigma_ij),j = 0  (LameEquation sets symmetry on)
        self.elasticity = LameEquation(domain)
        eoptions = self.elasticity.getSolverOptions()
        eoptions.setSolverMethod(elasticity_solver)
        # the elasticity operator is symmetric positive definite; for an iterative
        # method use an algebraic multigrid preconditioner (scales to fine meshes).
        # Tell AMG this is a `dim`-equations-per-node vector (elasticity) system so
        # it coarsens the block problem suitably rather than as a scalar PDE.
        if elasticity_solver != SolverOptions.DIRECT:
            eoptions.setPreconditioner(SolverOptions.AMG)
            eoptions.setTrilinosParameter("number of equations", dim)
            eoptions.setTrilinosParameter("verbosity", "none")
        self.helmholtz = None
        if self.localization_length is not None:
            c = self.localization_length ** 2 / 2.        # c = l**2/2
            self.helmholtz = LinearSinglePDE(domain, isComplex=False)
            self.helmholtz.setSymmetryOn()
            # the Helmholtz operator is symmetric positive definite, so PCG with
            # an algebraic multigrid preconditioner is used.
            options = self.helmholtz.getSolverOptions()
            options.setSolverMethod(SolverOptions.PCG)
            options.setPreconditioner(SolverOptions.AMG)
            # ebar - c*laplace(ebar) = etilde,  natural BC  grad(ebar).n = 0 (Eq. 20)
            self.helmholtz.setValue(A=c * kronecker(dim), D=1.)
        return self

    def getNonlocalStrain(self, etilde):
        """
        returns the non-local equivalent strain ebar as the solution of the
        implicit-gradient Helmholtz equation `ebar - c*laplace(ebar) = etilde`
        (Eq. 19) with natural boundary condition `grad(ebar).n = 0` (Eq. 20),
        interpolated to the element (`Function`) space of `etilde`.

        If no localization length was set, the local model is used and `etilde`
        is returned unchanged.
        """
        if self.helmholtz is None:
            return etilde
        self.helmholtz.setValue(Y=etilde)
        ebar = self.helmholtz.getSolution()
        return interpolate(ebar, etilde.getFunctionSpace())

    def getEquivalentStrain(self, eps):
        """
        modified von Mises equivalent strain (Eq. 9) for the symmetric strain
        tensor `eps` (rank-2 `Data`). Verified to return the axial strain for a
        uniaxial tension state.
        """
        g, nu = self.gamma, self.nu
        I1 = trace(eps)
        J2 = I1 ** 2 / 6. - inner(eps, eps) / 2.       # paper Eq. 10 convention
        c1 = (g - 1.) / (2. * g * (1. - 2. * nu))
        arg = ((g - 1.) ** 2 / (1. - 2. * nu) ** 2) * I1 ** 2 \
            - (12. * g / (1. + nu) ** 2) * J2
        return c1 * I1 + sqrt(clip(arg, minval=0.)) / (2. * g)

    def getDamage(self, kappa=None):
        """
        modified power-law damage D(kappa) (Eq. 8). Uses the stored history
        variable if `kappa` is not given. D=0 at kappa=kappa0 and D=1 for
        kappa>=kappa_c.
        """
        if kappa is None:
            kappa = self.kappa
        # `clip` takes scalar bounds only; `maximum`/`minimum` also accept the
        # `Data` bounds of a heterogeneous material.
        k = minimum(maximum(kappa, self.kappa0), self.kappa_c)
        return 1. - (self.kappa0 / k) ** self.beta \
            * ((self.kappa_c - k) / (self.kappa_c - self.kappa0)) ** self.alpha

    def update(self, eps_bar):
        """
        updates the history variable kappa from the (non-local) equivalent
        strain field `eps_bar` enforcing monotonicity (Eq. 7) and returns the
        updated damage D. On a single element pass the local equivalent strain
        `getEquivalentStrain(eps)`.
        """
        assert self.kappa is not None, "call initialize(domain) first."
        self.kappa = maximum(self.kappa, eps_bar)
        return self.getDamage()

    def getLameParameters(self, D):
        """
        returns the damaged isotropic Lame parameters (lambda_eff, mu_eff) for
        the damage field `D`. Since the damaged stiffness (1-D)*C0 is still
        isotropic, both parameters are the undamaged values scaled by (1-D),
        floored at `min_stiffness_ratio` to keep the equilibrium system
        non-singular as D -> 1. Ready to be used as the `lame_lambda`/`lame_mu`
        coefficients of a `LameEquation`.
        """
        g = clip(1. - D, minval=self.min_stiffness_ratio)
        return g * self.lam, g * self.mu

    def getStress(self, eps, D):
        """
        Cauchy stress sigma_ij = (1-D)*(lam*trace(eps)*delta_ij + 2*mu*eps_ij)
        for strain `eps` and damage `D`.
        """
        return (1. - D) * (self.lam * trace(eps) * kronecker(eps.getFunctionSpace())
                           + 2. * self.mu * eps)

    def getPorosity(self, eps, D=None):
        """
        auxiliary (one-way coupled) porosity field

            phi = phi0 + b(D) * trace(eps) + crack_porosity * D

        with the poroelastic volumetric term b(D)*trace(eps) (linear
        poroelasticity, dry case phi - phi0 = b*eps_v; Biot 1941, Coussy
        Poromechanics 2004) using a damage-degraded Biot coefficient

            b(D) = 1 - (1-D) K/Ks = biot + (1 - biot) * D

        (skeleton bulk modulus degraded by (1-D); biot = b0 = 1 - K/Ks is the
        intact value, b -> 1 as D -> 1). The optional crack_porosity*D term is a
        phenomenological dilatant crack porosity. NOTE b(D)*eps_v captures pore
        compaction/closure (in this damage-elasticity model eps_v is elastic, so
        under compression damage increases compaction); dilatant crack opening is
        represented by crack_porosity*D. Uses the stored damage if `D` is not
        given; clipped to [0, 1). Does NOT feed back into the deformation model.
        A permeability estimate follows e.g. Kozeny-Carman k = k0 phi^3/(1-phi)^2.
        """
        if D is None:
            D = self.D
        b = self.biot + (1. - self.biot) * D                 # b(D), b0 -> 1
        phi = self.porosity0 + b * trace(eps) + self.crack_porosity * D
        return clip(phi, minval=0., maxval=0.999)

    def solveLoadStep(self, tol=1e-6, max_iter=20, warn_on_nonconvergence=True):
        """
        advances the model by one displacement-controlled load step using the
        split-operator (staggered) scheme: the elasticity PDE `self.elasticity`,
        the implicit-gradient smoothing of the equivalent strain (if a
        localization length is set) and the damage update are solved alternately
        until the damage field converges.

        The equilibrium is solved incrementally: each staggered iteration solves
        for a displacement correction `du` from the internal-force residual of
        the stored total stress `self.stress`,

            -div(C(D):eps(du)) = div(self.stress),

        which enters the `LameEquation` as the prescribed stress `sigma =
        -self.stress` (the `X` coefficient). The load increment is applied as the
        constraint `du = r_total - self.u` on the constrained DOFs, where
        `r_total` is the (total) prescribed displacement set on the PDE. The
        correction is accumulated into the total displacement `self.u` and the
        stored stress is refreshed to `C(D):eps(u)` for the next residual, so the
        total secant equilibrium is preserved exactly as the damage `D` evolves.
        The traction-free faces are preserved by the natural boundary condition.

        The constraints (`q`, `r`) of `self.elasticity` must already be set for
        the current (total) load level. On entry the stored state `self.u`,
        `self.stress`, `self.D` and `self.kappa` is used as the starting point;
        on exit it is committed to the updated (non-decreasing damage) state.

        If the staggered loop does not reach `tol` within `max_iter` iterations
        the last iterate is committed anyway (the remaining imbalance is carried
        in `self.stress` and corrected at the next load step); the recommended
        remedy is to reduce the displacement-BC increment (sub-step the load
        level), which `runLoading` does automatically. The convergence flag is
        stored in `self.converged` and the final damage change in
        `self.last_change`.

        :param tol: convergence tolerance on the change of the damage field.
        :param max_iter: maximum number of staggered iterations.
        :param warn_on_nonconvergence: if True (default) a warning is logged when
                        the loop is exhausted without reaching `tol`. `runLoading`
                        sets this to False while sub-stepping so the automatic
                        bisection retries stay quiet.
        :return: (displacement `u`, strain `eps`, number of iterations used;
                 equal to `max_iter` when the loop did not converge).
        """
        assert self.elasticity is not None, "call initialize(domain) first."
        assert max_iter >= 1, "max_iter must be at least 1."
        D = self.D
        kappa_trial = self.kappa
        u = self.u
        stress = self.stress
        # total prescribed displacement for this load level (copied because the
        # PDE's r coefficient is overwritten with the increment below).
        r_total = self.elasticity.getCoefficient("r").copy()
        change = None
        for it in range(max_iter):
            lam_eff, mu_eff = self.getLameParameters(D)
            # incremental equilibrium: solve for the correction du.
            self.elasticity.setValue(lame_lambda=lam_eff, lame_mu=mu_eff,
                                     sigma=-stress, r=r_total - u)
            u = u + self.elasticity.getSolution()
            eps = symmetric(grad(u))
            etilde = self.getEquivalentStrain(eps)     # local equivalent strain
            ebar = self.getNonlocalStrain(etilde)      # implicit-gradient smoothing
            kappa_trial = maximum(self.kappa, ebar)
            D_trial = self.getDamage(kappa_trial)
            change = Lsup(D_trial - D)
            D = D_trial
            # refresh the stored total stress with the updated damage so it is the
            # residual base C(D):eps(u) of the next iteration / committed step.
            stress = self.getStress(eps, D)
            if change < tol:
                break
        self.last_change = change
        self.converged = bool(change < tol)
        if not self.converged and warn_on_nonconvergence:
            # max_iter exhausted without meeting tol: the committed (u, eps) are
            # slightly off-equilibrium for the final damage (the residual is
            # carried into self.stress and corrected at the next load step). The
            # better remedy is to reduce the displacement-BC increment (sub-step
            # this load level) so the staggered damage/equilibrium loop converges.
            self.logger.warning(
                "solveLoadStep did not converge: |dD|=%g >= tol=%g after %d "
                "iterations; consider reducing the load (BC) increment.",
                change, tol, max_iter)
        self.u = u                        # commit total displacement
        self.stress = stress              # commit total stress
        self.kappa = kappa_trial          # commit history (Eq. 7)
        self.D = D
        self.porosity = self.getPorosity(eps, D)   # auxiliary porosity
        return u, eps, it + 1

    def runLoading(self, set_bc, nsteps, callback=None, tol=1e-6, max_iter=20,
                   adaptive=True, min_increment=None, scale_to_onset=True):
        """
        runs a displacement-controlled loading path parameterised by a load
        factor in [0, 1], advancing the damage model with the split-operator
        scheme. The nominal increment is `1/nsteps`; the load level at factor
        `lambda` is prescribed by `set_bc(self.elasticity, lambda*nsteps, nsteps)`
        (so at the nominal factors `k/nsteps` this reduces to the plain
        `set_bc(pde, k, nsteps)`).

        For each (sub-)step:
          1. `set_bc` sets the constraints (`q`, `r`) of the elasticity PDE for
             the target (total) load level;
          2. `solveLoadStep` performs the staggered equilibrium/damage solve,
             updating the committed state (`self.u`, `self.stress`, `self.kappa`,
             `self.D`);
          3. the optional `callback(step, u, eps, self)` is invoked for output or
             recording of the solution `u`, strain `eps` and model state.

        Automatic sub-stepping (``adaptive=True``): if `solveLoadStep` does not
        converge, the committed state is rolled back, the increment is halved
        (bisection) and the (sub-)step is retried, down to `min_increment`. After
        a successful step the increment is grown back toward the nominal
        `1/nsteps`. If the (sub-)step still does not converge once the increment
        has been reduced to `min_increment`, the calculation is terminated by
        raising `RuntimeError`. With ``adaptive=False`` the classic fixed
        `nsteps` schedule is used (and `solveLoadStep` warns on non-convergence).

        First step scaled to damage onset (``scale_to_onset=True``): while the
        material is undamaged the response scales linearly with the load factor,
        and the (modified von Mises) equivalent strain is homogeneous of degree
        one in it. A single elastic solve at full load therefore gives the load
        factor `s = min(kappa0/ebar)` at which the non-local equivalent strain
        first reaches `kappa0` anywhere (damage onset; this is `kappa0/max(ebar)`
        for a uniform threshold). The committed state is set
        directly to this scaled elastic solution and stepping continues from `s`,
        skipping the purely-elastic sub-steps. This is only applied when starting
        from the undamaged state.

        Because `set_bc` is evaluated at intermediate load factors, it must
        accept a fractional `step` argument (the standard `u_bc*step/nsteps`
        form does).

        :param set_bc: callable `set_bc(pde, step, nsteps)` setting the load on
                       the elasticity PDE `pde` (= `self.elasticity`); may be
                       called with a fractional `step` when sub-stepping.
        :param nsteps: number of nominal load increments.
        :param callback: optional callable `callback(step, u, eps, model)`.
        :param tol: convergence tolerance passed to `solveLoadStep`.
        :param max_iter: maximum staggered iterations passed to `solveLoadStep`.
        :param adaptive: enable automatic bisection sub-stepping (default True).
        :param min_increment: smallest allowed load-factor increment (in [0, 1]);
                       the run is terminated with `RuntimeError` if a (sub-)step
                       at this increment does not converge. Defaults to
                       `1e-3 / nsteps`.
        :param scale_to_onset: if True (default) jump the first step elastically
                       to the damage-onset load factor (see above), applied only
                       when starting from the undamaged state.
        :return: list with the number of staggered iterations used per
                 (sub-)step (0 for the elastic onset step, if taken).
        :raises RuntimeError: if the minimum load increment is reached without
                       convergence.
        """
        assert self.elasticity is not None, "call initialize(domain) first."
        assert nsteps >= 1, "nsteps must be at least 1."
        iters = []
        nominal = 1. / nsteps
        if min_increment is None:
            min_increment = 1e-3 * nominal
        assert 0. < min_increment <= nominal, \
            "min_increment must be in (0, 1/nsteps]."
        dload = nominal                      # current load-factor increment
        load = 0.                            # committed load factor in [0, 1]
        step = 0
        # stop just short of 1 to avoid a spurious tiny final increment.
        load_end = 1. - 0.5 * min_increment

        # elastic jump to the damage-onset load factor (see docstring). Only
        # applied from the undamaged state; the elastic response scales linearly
        # with the load factor so a single full-load solve locates onset.
        if scale_to_onset and Lsup(self.D) == 0.:
            set_bc(self.elasticity, nsteps, nsteps)          # probe at factor 1
            r_full = self.elasticity.getCoefficient("r").copy()
            self.elasticity.setValue(lame_lambda=self.lam, lame_mu=self.mu,
                                     sigma=-self.stress, r=r_full)
            u_el = self.elasticity.getSolution()
            eps_el = symmetric(grad(u_el))
            ebar = self.getNonlocalStrain(self.getEquivalentStrain(eps_el))
            peak = sup(ebar)
            if peak > 0.:
                # smallest load factor at which kappa0 is reached anywhere. For a
                # uniform kappa0 this is just kappa0/sup(ebar); the pointwise form
                # also covers a heterogeneous threshold (layered specimen).
                s = inf(self.kappa0 / maximum(ebar, 1e-30 * peak))
                if s < 1.:
                    self.u = s * u_el
                    self.stress = self.getStress(s * eps_el, self.D)
                    load = s
                    step += 1
                    iters.append(0)
                    self.logger.info(
                        "scaled first step elastically to damage onset at load "
                        "factor %g.", s)
                    if callback is not None:
                        callback(step, self.u, s * eps_el, self)

        while load < load_end:
            target = min(load + dload, 1.)
            while True:
                if adaptive:
                    # snapshot the committed state so a non-converged attempt can
                    # be rolled back before retrying with a smaller increment.
                    saved = (self.u.copy(), self.stress.copy(),
                             self.kappa.copy(), self.D.copy())
                set_bc(self.elasticity, target * nsteps, nsteps)
                u, eps, nit = self.solveLoadStep(
                    tol=tol, max_iter=max_iter,
                    warn_on_nonconvergence=not adaptive)
                if self.converged or not adaptive:
                    break
                if dload <= min_increment:
                    # the minimum increment has been reached and still does not
                    # converge: terminate the calculation.
                    raise RuntimeError(
                        "load stepping failed to converge at load factor %g with "
                        "the minimum load increment %g (|dD|=%g >= tol=%g); "
                        "aborting." % (target, dload, self.last_change, tol))
                # reject: roll back and halve the increment (not below the
                # minimum), then retry.
                self.u, self.stress, self.kappa, self.D = saved
                dload = max(0.5 * dload, min_increment)
                target = min(load + dload, 1.)
                self.logger.info(
                    "load step did not converge; bisecting increment to %g "
                    "(target load factor %g).", dload, target)
            load = target
            step += 1
            iters.append(nit)
            if callback is not None:
                callback(step, u, eps, self)
            # grow the increment back toward the nominal size after success.
            if adaptive and dload < nominal:
                dload = min(2. * dload, nominal)
        return iters

    def dissipationRate(self, eps, D):
        """
        damage energy-release rate density Y = 1/2 eps:C0:eps (C0 the undamaged
        stiffness), the driving force of the dissipation. The incremental
        dissipation between two states is `integrate(Y * (D - D_old))`.
        """
        return 0.5 * (self.lam * trace(eps) ** 2 + 2. * self.mu * inner(eps, eps))

    def runLoadingDissipation(self, set_bc, callback=None, delta_tau=None,
                              tol=1e-6, max_iter=40, max_steps=300,
                              al_tol=1e-3, seed_factor=1.05, tau_growth=1.3,
                              min_delta_tau=None, stop_damage=0.999, stop=None):
        """
        dissipation-controlled (arc-length / path-following) load stepping that
        can traverse the post-peak softening branch on which displacement control
        (`runLoading`) stalls (snap-back).

        The load factor `lambda` becomes an unknown (the Dirichlet load is
        `r = lambda * r_ref`, `r_ref` the pattern set by `set_bc(pde, 1.0)`) and
        each step is controlled by prescribing the incremental dissipation
        `delta_tau = integrate(Y * (D - D_n))` (Y from `dissipationRate`). Since
        dissipation only increases (2nd law) it is a monotone path parameter even
        where `lambda` decreases, so the softening branch is followed.

        Each staggered corrector solves the residual correction `du_I` (r=0) and
        the load sensitivity `du_II` (r=r_ref, X=0) with the current secant
        stiffness, then a scalar solve gives `dlam` such that the total
        dissipation hits `delta_tau`; the damage is updated in the loop. The step
        is accepted when the equilibrium correction and the dissipation residual
        are below tolerance; otherwise `delta_tau` is halved and the step retried
        (down to `min_delta_tau`, below which `RuntimeError` is raised).

        Loading starts from the undamaged state: the first step scales elastically
        to just past damage onset (`seed_factor` times the onset factor) to seed a
        non-zero dissipation, after which `delta_tau` defaults to that seeding
        step's dissipation.

        :param set_bc: callable `set_bc(pde, factor)` setting `r = factor*r_ref`.
        :param callback: optional `callback(step, u, eps, model)`.
        :param delta_tau: prescribed dissipation increment per step (default: the
                        seeding step's dissipation).
        :param tol: equilibrium tolerance (relative Lsup of the correction).
        :param max_iter: max staggered correctors per step.
        :param max_steps: max number of (sub-)steps.
        :param al_tol: relative tolerance on the dissipation constraint.
        :param seed_factor: overshoot of the onset factor to initiate damage.
        :param tau_growth: growth factor applied to delta_tau after a good step.
        :param min_delta_tau: floor for delta_tau; RuntimeError below it.
        :param stop_damage: stop once the peak damage exceeds this (default
                        0.999); lower it to skip the stiff near-failure tail.
                        NOTE this tests the peak (`sup`) of the damage, so for a
                        statistically heterogeneous material a single failed weak
                        element ends the run; use `stop` instead there.
        :param stop: optional callable `stop(model)` evaluated before each step;
                        the run ends as soon as it returns True. Use it for
                        criteria that the peak damage cannot express, such as a
                        mean damage over a region or a drop off the peak load.
        :return: list of committed load factors `lambda` per step.
        """
        assert self.elasticity is not None, "call initialize(domain) first."
        assert Lsup(self.D) == 0., "runLoadingDissipation starts undamaged."
        set_bc(self.elasticity, 1.)                       # r_ref = r at factor 1
        r_ref = self.elasticity.getCoefficient("r").copy()
        zero_u = 0. * self.u
        zero_stress = 0. * self.stress

        # --- seed damage: elastic probe -> onset factor -> one step just past it
        self.elasticity.setValue(lame_lambda=self.lam, lame_mu=self.mu,
                                 sigma=-self.stress, r=r_ref)
        eps_el = symmetric(grad(self.elasticity.getSolution()))
        ebar_el = self.getNonlocalStrain(self.getEquivalentStrain(eps_el))
        peak = sup(ebar_el)
        assert peak > 0., "no straining under the reference load."
        # pointwise onset factor, see the note in `runLoading`.
        onset = inf(self.kappa0 / maximum(ebar_el, 1e-30 * peak))
        Dn = self.D.copy()                                # zero
        load = onset * seed_factor
        set_bc(self.elasticity, load)
        u, eps, _ = self.solveLoadStep(tol=tol, max_iter=max_iter,
                                       warn_on_nonconvergence=False)
        if delta_tau is None:
            delta_tau = integrate(self.dissipationRate(eps, self.D) * (self.D - Dn))
        if min_delta_tau is None:
            min_delta_tau = 1e-4 * delta_tau
        loads = [load]
        step = 1
        if callback is not None:
            callback(step, u, eps, self)
        dstep = load - onset                              # load-factor predictor

        # --- dissipation-controlled path following
        while step < max_steps and sup(self.D) <= stop_damage:
            if stop is not None and stop(self):
                self.logger.info("the `stop` criterion is met at step %d "
                                 "(load factor %g).", step, load)
                break
            saved = (self.u.copy(), self.stress.copy(),
                     self.kappa.copy(), self.D.copy())
            Dn = saved[3]
            kappa_n = saved[2]
            u, lam = self.u, load
            ok = False
            for it in range(max_iter):
                lam_eff, mu_eff = self.getLameParameters(self.D)
                self.elasticity.setValue(lame_lambda=lam_eff, lame_mu=mu_eff,
                                         sigma=-self.stress, r=zero_u)
                duI = self.elasticity.getSolution()       # residual correction
                self.elasticity.setValue(sigma=zero_stress, r=r_ref)
                duII = self.elasticity.getSolution()      # d u / d lambda

                def phi(dl):
                    e = symmetric(grad(u + duI + dl * duII))
                    eb = self.getNonlocalStrain(self.getEquivalentStrain(e))
                    Dt = self.getDamage(maximum(kappa_n, eb))
                    return integrate(self.dissipationRate(e, Dt)
                                     * (Dt - Dn)) - delta_tau

                dlam = _rootIncreasing(phi, max(abs(dstep), 1e-3 * abs(lam) + 1e-9),
                                       al_tol * delta_tau)
                u = u + duI + dlam * duII
                lam = lam + dlam
                eps = symmetric(grad(u))
                ebar = self.getNonlocalStrain(self.getEquivalentStrain(eps))
                self.kappa = maximum(kappa_n, ebar)
                self.D = self.getDamage(self.kappa)
                self.u = u
                self.stress = self.getStress(eps, self.D)
                eqn = Lsup(duI) / (Lsup(u) + 1e-30)
                gval = abs(integrate(self.dissipationRate(eps, self.D)
                                     * (self.D - Dn)) - delta_tau) / delta_tau
                if eqn < tol and gval < al_tol:
                    ok = True
                    break
            if not ok:
                # reject: roll back and reduce the dissipation increment
                self.u, self.stress, self.kappa, self.D = saved
                if delta_tau <= min_delta_tau:
                    raise RuntimeError(
                        "dissipation step failed to converge at the minimum "
                        "increment %g (load factor %g)." % (delta_tau, load))
                delta_tau *= 0.5
                self.logger.info("dissipation step did not converge; reducing "
                                 "delta_tau to %g.", delta_tau)
                continue
            dstep = lam - load
            load = lam
            step += 1
            loads.append(load)
            self.porosity = self.getPorosity(eps, self.D)   # auxiliary porosity
            if callback is not None:
                callback(step, self.u, eps, self)
            delta_tau *= tau_growth
        return loads

    @staticmethod
    def damageCurve(kappa, kappa0, kappa_c, alpha, beta):
        """
        `numpy` evaluation of the damage law D(kappa) (Eq. 8) for plotting /
        validation against Fig. 2. Accepts scalars or arrays.
        """
        k = np.clip(kappa, kappa0, kappa_c)
        return 1. - (kappa0 / k) ** beta \
            * ((kappa_c - k) / (kappa_c - kappa0)) ** alpha