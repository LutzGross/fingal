"""
Multi-element test of the smoothed continuum damage model of

    Mondal, Olsen-Kettle, Gross, Computers and Geotechnics 122 (2020) 103505.

A rectangular prism is loaded in displacement-controlled uniaxial compression
(smooth platens: only the axial displacement is prescribed on the top/bottom
faces, lateral expansion is free). A small weaker flaw (slightly reduced
stiffness E0, undamaged) is placed off-centre so that it strains more under load,
reaches the damage threshold first and nucleates damage that localises rather
than smearing uniformly -- this exercises the non-local implicit-gradient
smoothing, the incremental-residual equilibrium solve and the adaptive load
stepping that the single-element test (test1.py) cannot.

The mesh is resolved to two elements per localisation length l (the paper's
minimum to resolve the smoothing). The volume-averaged axial stress-strain
response is written to PNG, and the damage / displacement fields are written to
a silo file for visualisation (VisIt / ParaView).

Run:  PYTHONPATH=../../bin run-escript test2.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from esys.finley import Brick
from esys.escript import (Function, Solution, length, whereZero, whereNegative,
                          integrate, sup, Vector, Scalar)
from esys.escript.linearPDEs import SolverOptions
from esys.weipa import saveSilo

from fingal import SmoothDamageModel

# --------------------------------------------------------------------------
# material (as in test1.py: Table 1 damage law + Table 2 elastic properties)
# --------------------------------------------------------------------------
l_loc = 2.0e-3                                   # localisation length l = 2 mm
model = SmoothDamageModel(E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                          alpha=1.6, beta=0.01, gamma=841. / 250.,
                          localization_length=l_loc)

# --------------------------------------------------------------------------
# specimen and mesh: prism Lx x Ly x Lz, h = l/2 (two elements per l)
# --------------------------------------------------------------------------
Lx, Ly, Lz = 0.008, 0.008, 0.016                 # 8 x 8 x 16 mm
h = l_loc / 2.                                    # element size 1 mm
nx, ny, nz = int(round(Lx / h)), int(round(Ly / h)), int(round(Lz / h))
LOAD_DIR = 2                                      # compress along x3 (z)
domain = Brick(n0=nx, n1=ny, n2=nz, l0=Lx, l1=Ly, l2=Lz, order=1)
print(f"mesh: {nx} x {ny} x {nz} = {nx*ny*nz} hex elements (h = {h*1e3:g} mm, "
      f"{l_loc/h:g} elements per l)")

# iterative solvers for the multi-element problem (AMG multigrid preconditioner)
model.initialize(domain, elasticity_solver=SolverOptions.PCG)
eopts = model.elasticity.getSolverOptions()
eopts.setPreconditioner(SolverOptions.AMG)
eopts.setTolerance(1e-8)

# --------------------------------------------------------------------------
# weaker flaw: a small off-centre sphere with slightly reduced stiffness E0 (the
# Lame parameters lam, mu are scaled by `weak`; nu and the damage threshold
# kappa0 are unchanged). It is NOT pre-damaged -- under load the softer flaw
# strains more, reaches the equivalent-strain threshold kappa0 first and so
# nucleates damage there, which then localises.
# --------------------------------------------------------------------------
Xe = Function(domain).getX()                      # element-centre coordinates
flaw_centre = np.array([0.35 * Lx, 0.5 * Ly, 0.5 * Lz])
flaw_radius = 1.5 * h
weak = 0.8                                         # flaw stiffness = 0.8 * E0
inside = whereNegative(length(Xe - flaw_centre) - flaw_radius)
factor = 1. - (1. - weak) * inside                # `weak` inside, 1 outside
model.lam = model.lam * factor
model.mu = model.mu * factor
print(f"weaker flaw at {flaw_centre*1e3} mm, radius {flaw_radius*1e3:g} mm, "
      f"stiffness {weak:g}*E0 inside")

# --------------------------------------------------------------------------
# boundary conditions: uniaxial compression, smooth platens, minimal lateral
# pins to remove rigid-body motion (as in test1.py).
# --------------------------------------------------------------------------
x = domain.getX()
TOL = 1e-6 * h


def atpoint(p):
    return whereZero(length(x - np.array(p, dtype=float)), tol=TOL)


q = Vector(0., Solution(domain))
q[LOAD_DIR] = whereZero(x[2]) + whereZero(x[2] - Lz)
q[0] = atpoint((0., 0., 0.)) + atpoint((0., Ly, 0.))   # no x-transl / z-rot
q[1] = atpoint((0., 0., 0.)) + atpoint((Lx, 0., 0.))   # no y-transl
model.elasticity.setValue(q=q)

vol = integrate(Scalar(1., Function(domain)))

# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
# load up to the peak of the response; displacement-controlled loading deep into
# the post-peak softening branch is numerically stiff (snap-back) and makes the
# staggered scheme sub-step heavily -- arc-length control would be needed there.
u_bc_max = 1.25e-3 * Lz         # nominal global axial strain up to 1.25e-3
nsteps = 5
STAG_TOL, STAG_MAX = 1e-6, 30

hist = {"strain": [0.], "stress": [0.], "Dmax": [0.], "Dmean": [0.]}

# write the initial (pre-flaw) damage field as step 0
saveSilo("multi_element_step000", D=model.D)


def set_bc(pde, step, nsteps):
    r = Vector(0., Solution(domain))
    r[LOAD_DIR] = whereZero(x[2] - Lz) * (-u_bc_max * step / nsteps)
    pde.setValue(r=r)


def record(step, u, eps, model):
    D = model.D
    sig_zz = integrate(model.getStress(eps, D)[LOAD_DIR, LOAD_DIR]) / vol
    eps_zz = integrate(eps[LOAD_DIR, LOAD_DIR]) / vol
    hist["strain"].append(-eps_zz)
    hist["stress"].append(-sig_zz)
    hist["Dmax"].append(sup(D))
    hist["Dmean"].append(integrate(D) / vol)
    print(f"step {step:3d}: eps_zz={eps_zz: .3e} "
          f"sig_zz={sig_zz: .3e} Pa  D_max={sup(D): .4f}  "
          f"D_mean={integrate(D)/vol: .4f}", flush=True)
    # per-step field snapshot for visualising the localisation band
    saveSilo(f"multi_element_step{step:03d}", D=D, u=u,
             sig_zz=model.stress[LOAD_DIR, LOAD_DIR])


iters = model.runLoading(set_bc, nsteps, callback=record,
                         tol=STAG_TOL, max_iter=STAG_MAX)
print(f"executed {len(iters)} (sub-)steps; staggered iterations: {iters}")

# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------
plt.figure(figsize=(6, 4))
plt.plot(np.array(hist["strain"]) * 1e3, np.array(hist["stress"]) / 1e6, "b.-")
plt.xlabel("axial compressive strain  [1e-3]")
plt.ylabel("axial compressive stress [MPa]")
plt.title("Multi-element specimen: stress-strain")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("multi_element_response.png", dpi=120)
plt.close()
print("wrote multi_element_response.png")

# damage on the vertical mid-plane slice through the flaw (y = Ly/2), from the
# element-centre values -- shows the localisation band without a silo viewer.
coords = np.array(Function(domain).getX().toListOfTuples())
Dvals = np.array(model.D.toListOfTuples())
yc = coords[:, 1]
ytarget = yc[np.argmin(np.abs(yc - 0.5 * Ly))]     # nearest element-centre layer
mask = np.abs(yc - ytarget) < 0.25 * h
plt.figure(figsize=(4.2, 6))
# autoscale the colour range to the slice so the (mild, pre-peak) concentration
# at the weaker flaw is visible rather than washed out on a fixed [0, 1] scale.
sc = plt.scatter(coords[mask, 0] * 1e3, coords[mask, 2] * 1e3, c=Dvals[mask],
                 marker="s", s=260, cmap="inferno",
                 vmin=Dvals[mask].min(), vmax=Dvals[mask].max())
plt.colorbar(sc, label="damage  $D$")
plt.xlabel("x [mm]")
plt.ylabel("z [mm]")
plt.title(f"Damage on y = {ytarget*1e3:.1f} mm slice")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("multi_element_damage_slice.png", dpi=120)
plt.close()
print("wrote multi_element_damage_slice.png")

saveSilo("multi_element_damage", D=model.D, u=model.u,
         sig_zz=model.stress[LOAD_DIR, LOAD_DIR])
print("wrote multi_element_damage.silo")
print(f"peak stress = {max(hist['stress'])/1e6:.2f} MPa, "
      f"final D_max = {hist['Dmax'][-1]:.4f}, D_mean = {hist['Dmean'][-1]:.4f}")
