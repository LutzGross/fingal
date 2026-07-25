"""
Cubic specimen with a weak zone under uniaxial compression -- Section 3.1 /
Fig. 6-7 of Mondal, Olsen-Kettle, Gross, Computers and Geotechnics 122 (2020).

A 100 mm cube is compressed along x3 between CLAMPED (rough) platens
(u1=u2=u3=0 at x3=0; u1=u2=0, u3=-u_bc at x3=Lzh). A 5 mm cubic weak zone at the
centre of the front (x2-x3) face has a reduced stiffness (E=18 GPa, nu=0.26 vs
E=36 GPa, nu=0.18 outside). Damage nucleates at the weak zone and propagates as
conjugate shear bands (the X-pattern of Fig. 7b).

This is a COARSE, qualitative reproduction: the paper uses ~40 million elements
(h=1 mm) on HPC; here the mesh is set by the CLI argument. The Weibull
heterogeneity of the outer Young's modulus is omitted (its fine detail is smeared
by the non-local length l=2 mm; the weak zone is the dominant trigger).

Run:  PYTHONPATH=../../bin run-escript test5.py [n_per_side] [stop_damage]
"""
import os
import sys
import numpy as np

from esys.finley import Brick
from esys.escript import (Function, Solution, whereZero, whereNegative,
                          wherePositive, integrate, sup, Vector, Scalar, kronecker,
                          trace)
from esys.escript.linearPDEs import SolverOptions
from esys.weipa import saveSilo

from fingal import SmoothDamageModel

n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
stop_damage = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
RESULTS = f"cube_n{n}"
os.makedirs(RESULTS, exist_ok=True)

L = 0.1                                        # 100 mm cube
LOAD_DIR = 2
l_loc = 2.0e-3
model = SmoothDamageModel(E0=36.e9, nu=0.18, kappa0=2.0e-4, kappa_c=1.0e-3,
                          alpha=1.6, beta=0.01, gamma=841. / 250.,
                          localization_length=l_loc,
                          porosity0=0.1, crack_porosity=0.1)  # illustrative
domain = Brick(n0=n, n1=n, n2=n, l0=L, l1=L, l2=L, order=1)
h = L / n
print(f"mesh: {n}^3 = {n**3} elem, h={h*1e3:g} mm ({l_loc/h:g} elem/l)", flush=True)
model.initialize(domain, elasticity_solver=SolverOptions.PCG)
model.elasticity.getSolverOptions().setTolerance(1e-8)

# --- weak zone: 5 mm cube at the centre of the front (x1=0) x2-x3 face -------
# reduced Young's modulus and raised Poisson ratio inside; set the (spatial)
# Lame parameters directly from the E(x), nu(x) fields.
Xe = Function(domain).getX()
inside = (whereNegative(Xe[0] - 5e-3)
          * whereNegative((Xe[1] - L / 2) ** 2 - 2.5e-3 ** 2)
          * whereNegative((Xe[2] - L / 2) ** 2 - 2.5e-3 ** 2))
E = 36.e9 - (36.e9 - 18.e9) * inside
nu = 0.18 + (0.26 - 0.18) * inside
model.lam = E * nu / ((1. + nu) * (1. - 2. * nu))
model.mu = E / (2. * (1. + nu))
print(f"weak zone volume fraction = {float(integrate(inside)/L**3):.2e}", flush=True)

# --- clamped (rough) platens, uniaxial compression along x3 ------------------
x = domain.getX()
bottom, top = whereZero(x[2]), whereZero(x[2] - L)
q = Vector(0., Solution(domain))
q[0] = bottom + top
q[1] = bottom + top
q[2] = bottom + top
model.elasticity.setValue(q=q)
u_bc_ref = 4.0e-3 * L           # reference top displacement (load factor 1)
vol = integrate(Scalar(1., Function(domain)))


def set_bc(pde, factor):
    r = Vector(0., Solution(domain))
    r[2] = top * (-u_bc_ref * factor)
    pde.setValue(r=r)


# --- stream stress-strain history + per-step damage field -------------------
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
    print(f"step {step:3d}: eps={strain: .3e} sig={stress/1e6: .2f}MPa "
          f"Dmax={dmax: .4f} Dmean={dmean: .4f}", flush=True)
    saveSilo(os.path.join(RESULTS, f"step{step:03d}"), D=D, u=u,
             sig_zz=model.stress[LOAD_DIR, LOAD_DIR], phi=model.porosity)


loads = model.runLoadingDissipation(set_bc, callback=record, tol=1e-6,
                                    max_iter=40, max_steps=120, al_tol=1e-3,
                                    stop_damage=stop_damage)
hist_csv.close()
print(f"done: {len(loads)} steps, peak {max(hist['stress'])/1e6:.2f} MPa, "
      f"Dmax {hist['Dmax'][-1]:.3f}", flush=True)
print(f"visualise {RESULTS}/step*.silo in VisIt; plot with: python3 plot.py {RESULTS}",
      flush=True)
