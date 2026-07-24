"""
Dissipation-controlled (arc-length) loading of the smoothed continuum damage
model -- the post-peak counterpart of test2.py.

The same weaker-flaw prism as test2.py is driven with
SmoothDamageModel.runLoadingDissipation, which makes the load factor an unknown
and prescribes the incremental dissipation per step. Because dissipation only
increases (2nd law) it is a monotone path parameter even where the load factor
decreases, so the softening branch -- on which displacement control (test2.py)
stalls with heavy sub-stepping -- can be traversed and a localisation band
develops.

Run:  PYTHONPATH=../../bin run-escript test3.py [elements_per_l] [stop_damage]

  elements_per_l : mesh resolution as elements per localisation length l
                   (default 2). Increasing it (e.g. 3, 4) refines the mesh for a
                   mesh-objectivity study -- the peak stress and band width
                   should stay ~constant because l fixes the band width.
  stop_damage    : stop once peak damage exceeds this (default 0.999); lower it
                   (e.g. 0.9) on fine meshes to skip the stiff near-failure tail.
"""
import sys
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

n_per_l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
stop_damage = float(sys.argv[2]) if len(sys.argv) > 2 else 0.999
tag = "" if n_per_l == 2 else f"_{n_per_l}pl"       # output-name suffix

l_loc = 2.0e-3
model = SmoothDamageModel(E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                          alpha=1.6, beta=0.01, gamma=841. / 250.,
                          localization_length=l_loc)

Lx, Ly, Lz = 0.008, 0.008, 0.016
h = l_loc / n_per_l
nx, ny, nz = int(round(Lx / h)), int(round(Ly / h)), int(round(Lz / h))
LOAD_DIR = 2
domain = Brick(n0=nx, n1=ny, n2=nz, l0=Lx, l1=Ly, l2=Lz, order=1)
print(f"mesh: {nx} x {ny} x {nz} = {nx*ny*nz} hex elements "
      f"(h = {h*1e3:g} mm, {n_per_l} elements per l)", flush=True)

# PCG + AMG multigrid (set up by initialize); tighten the tolerance
model.initialize(domain, elasticity_solver=SolverOptions.PCG)
model.elasticity.getSolverOptions().setTolerance(1e-9)

# weaker flaw (reduced E0), as in test2.py. The flaw is a fixed PHYSICAL size
# (independent of h) so the mesh-refinement study is meaningful.
Xe = Function(domain).getX()
flaw_centre = np.array([0.35 * Lx, 0.5 * Ly, 0.5 * Lz])
flaw_radius = 0.75 * l_loc                          # = 1.5 * h at 2 elements per l
weak = 0.7                                          # flaw stiffness = 0.7 * E0
inside = whereNegative(length(Xe - flaw_centre) - flaw_radius)
factor = 1. - (1. - weak) * inside
model.lam = model.lam * factor
model.mu = model.mu * factor
print(f"weaker flaw, stiffness {weak:g}*E0", flush=True)

x = domain.getX()
TOL = 1e-6 * h
atpoint = lambda p: whereZero(length(x - np.array(p, float)), tol=TOL)
q = Vector(0., Solution(domain))
q[LOAD_DIR] = whereZero(x[2]) + whereZero(x[2] - Lz)
q[0] = atpoint((0., 0., 0.)) + atpoint((0., Ly, 0.))
q[1] = atpoint((0., 0., 0.)) + atpoint((Lx, 0., 0.))
model.elasticity.setValue(q=q)

vol = integrate(Scalar(1., Function(domain)))
u_bc_ref = 4.0e-3 * Lz          # reference top displacement (load factor 1)


def set_bc(pde, factor):
    r = Vector(0., Solution(domain))
    r[LOAD_DIR] = whereZero(x[2] - Lz) * (-u_bc_ref * factor)
    pde.setValue(r=r)


hist = {"strain": [], "stress": [], "Dmax": [], "Dmean": []}


def record(step, u, eps, model):
    D = model.D
    hist["strain"].append(-integrate(eps[LOAD_DIR, LOAD_DIR]) / vol)
    hist["stress"].append(-integrate(model.getStress(eps, D)[LOAD_DIR, LOAD_DIR]) / vol)
    hist["Dmax"].append(sup(D))
    hist["Dmean"].append(integrate(D) / vol)
    print(f"step {step:3d}: eps_zz={-hist['strain'][-1]: .3e} "
          f"sig_zz={-hist['stress'][-1]: .3e} Pa  D_max={sup(D): .4f}  "
          f"D_mean={hist['Dmean'][-1]: .4f}", flush=True)
    saveSilo(f"m3{tag}_step{step:03d}", D=D, u=u,
             sig_zz=model.stress[LOAD_DIR, LOAD_DIR])


loads = model.runLoadingDissipation(set_bc, callback=record, tol=1e-6,
                                    max_iter=40, max_steps=80, al_tol=1e-3,
                                    stop_damage=stop_damage)
print(f"executed {len(loads)} steps; load factor {loads[0]:.3f} -> {loads[-1]:.3f} "
      f"(min {min(loads):.3f})", flush=True)

# stress-strain (should now include the descending softening branch)
plt.figure(figsize=(6, 4))
plt.plot(np.array(hist["strain"]) * 1e3, np.array(hist["stress"]) / 1e6, "b.-")
plt.xlabel("axial compressive strain  [1e-3]")
plt.ylabel("axial compressive stress [MPa]")
plt.title(f"Dissipation-controlled stress-strain ({n_per_l} elem/l)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"m3_response{tag}.png", dpi=120)
plt.close()
print(f"wrote m3_response{tag}.png", flush=True)

# damage mid-plane slice
coords = np.array(Function(domain).getX().toListOfTuples())
Dvals = np.array(model.D.toListOfTuples())
yc = coords[:, 1]
ytarget = yc[np.argmin(np.abs(yc - 0.5 * Ly))]
mask = np.abs(yc - ytarget) < 0.25 * h
plt.figure(figsize=(4.2, 6))
sc = plt.scatter(coords[mask, 0] * 1e3, coords[mask, 2] * 1e3, c=Dvals[mask],
                 marker="s", s=260 * (2. / n_per_l) ** 2, cmap="inferno",
                 vmin=0., vmax=1.)
plt.colorbar(sc, label="damage  $D$")
plt.xlabel("x [mm]")
plt.ylabel("z [mm]")
plt.title(f"Damage slice, {n_per_l} elem/l  (y = {ytarget*1e3:.1f} mm)")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig(f"m3_damage_slice{tag}.png", dpi=120)
plt.close()
print(f"wrote m3_damage_slice{tag}.png", flush=True)
print(f"peak stress = {max(hist['stress'])/1e6:.2f} MPa, "
      f"final D_max = {hist['Dmax'][-1]:.4f}", flush=True)
