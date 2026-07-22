"""
Single finite-element test of the smoothed continuum damage model of

    Mondal, Olsen-Kettle, Gross, Computers and Geotechnics 122 (2020) 103505.

Two checks are performed:

  (1) Constitutive check: the damage law D(kappa) (Eq. 8) is evaluated over
      [kappa0, kappa_c] and compared to Fig. 2 of the paper.

  (2) FEM check: a single tri-linear hexahedral element is loaded in
      displacement-controlled uniaxial compression with free lateral faces
      (minimal corner pins allow Poisson expansion). A staggered loop updates
      damage and re-solves equilibrium, producing the axial stress-strain and
      damage-strain response.

On a single element the non-local Helmholtz smoothing is the identity
(uniform field, Neumann boundary), so this exercises the full local model and
the damage/elasticity coupling before moving to real meshes.

Run:  PYTHONPATH=../../bin run-escript test1.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from esys.finley import Brick
from esys.escript import (Function, Solution, length, whereZero,
                          integrate, Vector, Scalar)
from esys.escript.linearPDEs import SolverOptions

from fingal import SmoothDamageModel

# --------------------------------------------------------------------------
# material (Table 1: damage law; Table 2: elastic properties, strength ratio)
# --------------------------------------------------------------------------
model = SmoothDamageModel(E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                          alpha=1.6, beta=0.01, gamma=841. / 250.,
                          localization_length=2.0e-3)   # l = 2 mm (Table 2)

# ==========================================================================
# (1) constitutive check -- reproduce Fig. 2 (damage vs equivalent strain)
# ==========================================================================
kappa = np.linspace(0., 1.2 * model.kappa_c, 601)
Dcurve = model.damageCurve(kappa, model.kappa0, model.kappa_c,
                           model.alpha, model.beta)
plt.figure(figsize=(6, 4))
plt.plot(kappa, Dcurve, "r-", lw=2)
plt.xlabel("equivalent strain  $\\tilde{\\varepsilon}=\\kappa$")
plt.ylabel("damage  $D$")
plt.title("Damage versus strain (cf. Fig. 2)")
plt.xlim(0., 1.2 * model.kappa_c)
plt.ylim(0., 1.05)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("damage_vs_strain.png", dpi=120)
plt.close()
print("wrote damage_vs_strain.png")

# ==========================================================================
# (2) single-element FEM uniaxial compression
# ==========================================================================
L = 0.005                       # element edge length [m] (5 mm)
LOAD_DIR = 2                     # compress along x3 (z)
domain = Brick(n0=1, n1=1, n2=1, l0=L, l1=L, l2=L, order=1)
# the model owns the elasticity LameEquation (self.elasticity); use a direct
# solver for this tiny single-element problem.
model.initialize(domain, elasticity_solver=SolverOptions.DIRECT)

# --- boundary conditions: uniaxial compression, free lateral faces --------
x = domain.getX()
TOL = 1e-8 * L


def atpoint(p):
    """nodal mask at corner p=(x,y,z)."""
    return whereZero(length(x - np.array(p, dtype=float)), tol=TOL)


# constrain u_z on bottom (z=0 -> 0) and top (z=L -> -u_bc); pin x,y minimally
q = Vector(0., Solution(domain))
q[LOAD_DIR] = whereZero(x[2]) + whereZero(x[2] - L)
q[0] = atpoint((0., 0., 0.)) + atpoint((0., L, 0.))   # no x-transl / z-rot
q[1] = atpoint((0., 0., 0.)) + atpoint((L, 0., 0.))   # no y-transl
model.elasticity.setValue(q=q)

vol = integrate(Scalar(1., Function(domain)))         # element volume

# --- displacement-controlled loading loop ---------------------------------
u_bc_max = 4.0e-3 * L           # axial strain up to ~gamma*kappa_c (full failure)
nsteps = 80
STAG_TOL = 1e-6
STAG_MAX = 20

hist = {"strain": [], "stress": [], "damage": [], "eqstrain": []}


def set_bc(pde, step, nsteps):
    """apply the axial compression displacement for the given load step."""
    u_bc = u_bc_max * step / nsteps
    r = Vector(0., Solution(domain))
    r[LOAD_DIR] = whereZero(x[2] - L) * (-u_bc)
    pde.setValue(r=r)


def record(step, u, eps, model):
    """volume-average and store the single-element response of a load step."""
    D = model.D
    eq = model.getEquivalentStrain(eps)
    sig = model.getStress(eps, D)
    sig_zz = integrate(sig[LOAD_DIR, LOAD_DIR]) / vol
    eps_zz = integrate(eps[LOAD_DIR, LOAD_DIR]) / vol
    eq_avg = integrate(eq) / vol
    D_avg = integrate(D) / vol
    hist["strain"].append(-eps_zz)          # compressive strain (positive)
    hist["stress"].append(-sig_zz)          # compressive stress (positive)
    hist["damage"].append(D_avg)
    hist["eqstrain"].append(eq_avg)
    print(f"step {step:3d}: eps_zz={eps_zz: .3e} "
          f"eq={eq_avg: .3e} D={D_avg: .4f} sig_zz={sig_zz: .3e} Pa")


model.runLoading(set_bc, nsteps, callback=record,
                 tol=STAG_TOL, max_iter=STAG_MAX)

# --- plot single-element response -----------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(np.array(hist["strain"]), np.array(hist["stress"]) / 1e6, "b.-")
ax[0].set_xlabel("axial compressive strain")
ax[0].set_ylabel("axial compressive stress [MPa]")
ax[0].set_title("Single-element stress-strain")
ax[0].grid(True, alpha=0.3)

ax[1].plot(hist["eqstrain"], hist["damage"], "r.-")
ax[1].set_xlabel("equivalent strain  $\\tilde{\\varepsilon}$")
ax[1].set_ylabel("damage  $D$")
ax[1].set_title("Single-element damage evolution")
ax[1].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("single_element_response.png", dpi=120)
plt.close()
print("wrote single_element_response.png")
print(f"peak stress = {max(hist['stress']) / 1e6:.2f} MPa, "
      f"final damage = {hist['damage'][-1]:.4f}")