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
import os
import sys
import numpy as np

from esys.finley import Brick
from esys.escript import (Function, Solution, length, whereZero, whereNegative,
                          integrate, sup, Vector, Scalar)
from esys.escript.linearPDEs import SolverOptions
from esys.weipa import saveSilo

from fingal import SmoothDamageModel

n_per_l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
stop_damage = float(sys.argv[2]) if len(sys.argv) > 2 else 0.999
# all outputs of this case go into their own directory
RESULTS = f"dissipation_{n_per_l}pl"
os.makedirs(RESULTS, exist_ok=True)

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


# fixed mid-plane (y = Ly/2) slice geometry for the streamed slice.csv
_coords = np.array(Function(domain).getX().toListOfTuples())
_yc = _coords[:, 1]
_yt = _yc[np.argmin(np.abs(_yc - 0.5 * Ly))]
_smask = np.abs(_yc - _yt) < 0.25 * h
_sx, _sz = _coords[_smask, 0] * 1e3, _coords[_smask, 2] * 1e3

# stream the plot data to CSV as the run proceeds, so a plot can be made (with
# plot.py) even if the run is terminated early.
hist = {"stress": [], "Dmax": []}
hist_csv = open(os.path.join(RESULTS, "history.csv"), "w")
hist_csv.write("step,strain,stress,Dmax,Dmean\n")
hist_csv.flush()


def record(step, u, eps, model):
    D = model.D
    strain = -float(integrate(eps[LOAD_DIR, LOAD_DIR]) / vol)
    stress = -float(integrate(model.getStress(eps, D)[LOAD_DIR, LOAD_DIR]) / vol)
    dmax, dmean = sup(D), float(integrate(D) / vol)
    hist["stress"].append(stress)
    hist["Dmax"].append(dmax)
    hist_csv.write(f"{step},{strain:.8e},{stress:.8e},{dmax:.6f},{dmean:.6f}\n")
    hist_csv.flush()
    print(f"step {step:3d}: eps_zz={-strain: .3e} sig_zz={-stress: .3e} Pa  "
          f"D_max={dmax: .4f}  D_mean={dmean: .4f}", flush=True)
    saveSilo(os.path.join(RESULTS, f"step{step:03d}"), D=D, u=u,
             sig_zz=model.stress[LOAD_DIR, LOAD_DIR])
    # current mid-plane damage slice (overwritten so the latest always survives)
    Dv = np.array(D.toListOfTuples())[_smask]
    with open(os.path.join(RESULTS, "slice.csv"), "w") as f:
        f.write("x_mm,z_mm,D\n")
        for xv, zv, dv in zip(_sx, _sz, Dv):
            f.write(f"{xv:.4f},{zv:.4f},{dv:.6f}\n")


loads = model.runLoadingDissipation(set_bc, callback=record, tol=1e-6,
                                    max_iter=40, max_steps=80, al_tol=1e-3,
                                    stop_damage=stop_damage)
hist_csv.close()
print(f"executed {len(loads)} steps; load factor {loads[0]:.3f} -> {loads[-1]:.3f} "
      f"(min {min(loads):.3f})", flush=True)
print(f"peak stress = {max(hist['stress'])/1e6:.2f} MPa, "
      f"final D_max = {hist['Dmax'][-1]:.4f}", flush=True)
print(f"plot with:  python3 plot.py {RESULTS}", flush=True)
