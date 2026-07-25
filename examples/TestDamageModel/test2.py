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
import os
import numpy as np

from esys.finley import Brick
from esys.escript import (Function, Solution, length, whereZero, whereNegative,
                          integrate, sup, Vector, Scalar)
from esys.escript.linearPDEs import SolverOptions
from esys.weipa import saveSilo

from fingal import SmoothDamageModel

RESULTS = "multi_element"                 # all outputs of this case go here
os.makedirs(RESULTS, exist_ok=True)

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

# fixed mid-plane (y = Ly/2) slice geometry for the streamed slice.csv
_coords = np.array(Function(domain).getX().toListOfTuples())
_yc = _coords[:, 1]
_yt = _yc[np.argmin(np.abs(_yc - 0.5 * Ly))]
_smask = np.abs(_yc - _yt) < 0.25 * h
_sx, _sz = _coords[_smask, 0] * 1e3, _coords[_smask, 2] * 1e3


def write_slice(D):
    """overwrite slice.csv with the current mid-plane damage (survives abort)."""
    Dv = np.array(D.toListOfTuples())[_smask]
    with open(os.path.join(RESULTS, "slice.csv"), "w") as f:
        f.write("x_mm,z_mm,D\n")
        for xv, zv, dv in zip(_sx, _sz, Dv):
            f.write(f"{xv:.4f},{zv:.4f},{dv:.6f}\n")


# stream the plot data to CSV as the run proceeds (plot with plot.py); seed the
# unloaded origin so the elastic branch is drawn from (0, 0).
hist = {"stress": [0.], "Dmax": [0.]}
hist_csv = open(os.path.join(RESULTS, "history.csv"), "w")
hist_csv.write("step,strain,stress,Dmax,Dmean\n")
hist_csv.write("0,0.0,0.0,0.0,0.0\n")
hist_csv.flush()
saveSilo(os.path.join(RESULTS, "step000"), D=model.D)
write_slice(model.D)


def set_bc(pde, step, nsteps):
    r = Vector(0., Solution(domain))
    r[LOAD_DIR] = whereZero(x[2] - Lz) * (-u_bc_max * step / nsteps)
    pde.setValue(r=r)


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
    write_slice(D)


iters = model.runLoading(set_bc, nsteps, callback=record,
                         tol=STAG_TOL, max_iter=STAG_MAX)
hist_csv.close()
print(f"executed {len(iters)} (sub-)steps; staggered iterations: {iters}")
saveSilo(os.path.join(RESULTS, "multi_element_damage"), D=model.D, u=model.u,
         sig_zz=model.stress[LOAD_DIR, LOAD_DIR])
print(f"peak stress = {max(hist['stress'])/1e6:.2f} MPa, "
      f"final D_max = {hist['Dmax'][-1]:.4f}", flush=True)
print(f"plot with:  python3 plot.py {RESULTS}", flush=True)
