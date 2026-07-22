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

import numpy as np
from esys.escript import (Function, Scalar, kronecker, trace, inner,
                          sqrt, clip, maximum, symmetric, grad, Lsup, interpolate)
from esys.escript.linearPDEs import LinearSinglePDE, LameEquation, SolverOptions

__all__ = ['SmoothDamageModel']


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
    """
    def __init__(self, E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                 alpha=1.6, beta=0.01, gamma=841. / 250., sigma_t=None,
                 min_stiffness_ratio=1.e-3, localization_length=None):
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
        """
        assert 0. <= nu < 0.5, "Poisson ratio must be in [0, 0.5)."
        assert kappa_c > kappa0 > 0., "need kappa_c > kappa0 > 0."
        assert gamma >= 1., "compressive strength must not be below tensile strength."
        self.E0 = E0
        self.nu = nu
        if sigma_t is not None:
            kappa0 = sigma_t / E0
        self.kappa0 = kappa0
        self.kappa_c = kappa_c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.min_stiffness_ratio = min_stiffness_ratio
        self.localization_length = localization_length
        # isotropic Lame parameters of the undamaged material; the full 4th
        # order stiffness tensor is `isotropicStiffnessTensor(E0, nu)`.
        self.lam = E0 * nu / ((1. + nu) * (1. - 2. * nu))
        self.mu = E0 / (2. * (1. + nu))
        self.kappa = None
        self.elasticity = None
        self.helmholtz = None

    def initialize(self, domain, elasticity_solver=SolverOptions.PCG):
        """
        prepares the model on `domain`: sets the history variable kappa to the
        threshold value kappa0 and the damage D to zero on the element
        (`Function`) function space, builds the `LameEquation` for the
        elasticity problem (available as `self.elasticity`) and, if a
        localization length is set, the Helmholtz PDE for the non-local
        smoothing. Returns self.

        The caller sets the constraints (`q`, `r`) on `self.elasticity` for the
        loading; the damaged Lame parameters are (re)set by `solveLoadStep`.

        :param elasticity_solver: solver method for the elasticity PDE.
        """
        self.kappa = Scalar(self.kappa0, Function(domain))
        self.D = Scalar(0., Function(domain))
        # elasticity problem  -(sigma_ij),j = 0  (LameEquation sets symmetry on)
        self.elasticity = LameEquation(domain)
        self.elasticity.getSolverOptions().setSolverMethod(elasticity_solver)
        self.helmholtz = None
        if self.localization_length is not None:
            dim = domain.getDim()
            c = self.localization_length ** 2 / 2.        # c = l**2/2
            self.helmholtz = LinearSinglePDE(domain, isComplex=False)
            self.helmholtz.setSymmetryOn()
            self.helmholtz.getSolverOptions().setSolverMethod(SolverOptions.PCG)
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
        k = clip(kappa, minval=self.kappa0, maxval=self.kappa_c)
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
        return (1. - D) * (self.lam * trace(eps) * kronecker(3) + 2. * self.mu * eps)

    def solveLoadStep(self, tol=1e-6, max_iter=20):
        """
        advances the model by one displacement-controlled load step using the
        split-operator (staggered) scheme: the elasticity PDE `self.elasticity`,
        the implicit-gradient smoothing of the equivalent strain (if a
        localization length is set) and the damage update are solved alternately
        until the damage field converges.

        The constraints (`q`, `r`) of `self.elasticity` must already be set for
        the current load level; its Lame coefficients are (re)set here from the
        current damage. On entry the stored damage `self.D` and history
        `self.kappa` are used as the starting point; on exit the history is
        committed and `self.D` holds the (non-decreasing) updated damage.

        :param tol: convergence tolerance on the change of the damage field.
        :param max_iter: maximum number of staggered iterations.
        :return: (displacement `u`, strain `eps`, number of iterations used).
        """
        assert self.elasticity is not None, "call initialize(domain) first."
        D = self.D
        kappa_trial = self.kappa
        for it in range(max_iter):
            lam_eff, mu_eff = self.getLameParameters(D)
            self.elasticity.setValue(lame_lambda=lam_eff, lame_mu=mu_eff)
            u = self.elasticity.getSolution()
            eps = symmetric(grad(u))
            etilde = self.getEquivalentStrain(eps)     # local equivalent strain
            ebar = self.getNonlocalStrain(etilde)      # implicit-gradient smoothing
            kappa_trial = maximum(self.kappa, ebar)
            D_trial = self.getDamage(kappa_trial)
            change = Lsup(D_trial - D)
            D = D_trial
            if change < tol:
                break
        self.kappa = kappa_trial          # commit history (Eq. 7)
        self.D = D
        return u, eps, it + 1

    def runLoading(self, set_bc, nsteps, callback=None, tol=1e-6, max_iter=20):
        """
        runs a displacement-controlled loading path of `nsteps` increments,
        advancing the damage model with the split-operator scheme at each step.

        For each step k = 1 .. nsteps:
          1. `set_bc(self.elasticity, k, nsteps)` sets the constraints (`q`,
             `r`) of the elasticity PDE for the current load level;
          2. `solveLoadStep` performs the staggered equilibrium/damage solve,
             updating `self.kappa` and `self.D`;
          3. the optional `callback(k, u, eps, self)` is invoked for output or
             recording of the solution `u`, strain `eps` and model state.

        :param set_bc: callable `set_bc(pde, step, nsteps)` setting the load on
                       the elasticity PDE `pde` (= `self.elasticity`).
        :param nsteps: number of load increments.
        :param callback: optional callable `callback(step, u, eps, model)`.
        :param tol: convergence tolerance passed to `solveLoadStep`.
        :param max_iter: maximum staggered iterations passed to `solveLoadStep`.
        :return: list with the number of staggered iterations used per step.
        """
        assert self.elasticity is not None, "call initialize(domain) first."
        iters = []
        for step in range(1, nsteps + 1):
            set_bc(self.elasticity, step, nsteps)
            u, eps, nit = self.solveLoadStep(tol=tol, max_iter=max_iter)
            iters.append(nit)
            if callback is not None:
                callback(step, u, eps, self)
        return iters

    @staticmethod
    def damageCurve(kappa, kappa0, kappa_c, alpha, beta):
        """
        `numpy` evaluation of the damage law D(kappa) (Eq. 8) for plotting /
        validation against Fig. 2. Accepts scalars or arrays.
        """
        k = np.clip(kappa, kappa0, kappa_c)
        return 1. - (kappa0 / k) ** beta \
            * ((kappa_c - k) / (kappa_c - kappa0)) ** alpha