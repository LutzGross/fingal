#!/usr/bin/python3
from esys.escript import *
from esys.finley import ReadMesh
from esys.weipa import saveSilo
from fingal import setupERTPDE, SIPSolver
from esys.escript.linearPDEs import SolverOptions
import time, json

import argparse


parser = argparse.ArgumentParser(description='run SIP solver test for unstructured mesh.')
parser.add_argument(dest='meshfile', metavar='meshfile', type=str, help='fly file')
#parser.add_argument('--silo', '-s',  dest='silo', metavar='SILO', default=None, help="silo file to write.")
parser.add_argument('--contrast', '-C',  dest='contrast', metavar='CONTRAST', type=float, default=1, help="raise factor for anomaly conductivity vs. background.")
parser.add_argument('--ratio', '-R',  dest='ratio', metavar='RATIO', type=float, default=0, help="ratio imaginary over real part of conductivity in anomaly.")
parser.add_argument('--silent', '-S',  dest='silent', action='store_true', default=False, help="no info print")
parser.add_argument('--test', '-t',  dest='test', action='store_true', default=False, help="test against a direct solver solution")

args = parser.parse_args()

#-------------------------------------------------------------
SIGMA_REAL = 1.
SIGMA_RATIO = 0.01 # ratio real to imag
ANOMALY_SIGMA_REAL = SIGMA_REAL * args.contrast
ANOMALY_SIGMA_RATIO = args.ratio # ratio real to imag
#--------------------------------------------------------------

SRC_LOC = (0.,0.,0.)
SRC_I = 1.

dom = ReadMesh(args.meshfile)
print("mesh read from ", args.meshfile)

# ... mask of the insulated surface.
surface_mask = Scalar(0., FunctionOnBoundary(dom))
surface_mask.setTaggedValue("top", 1)
# ... set-up conductivity distribution:
# ...... the overall conductiviry also acts as primary conductivity
sigma_bg = SIGMA_REAL * (1 + 1j * SIGMA_RATIO)
# ....... conductivity in the anomaly
sigma_anomaly = ANOMALY_SIGMA_REAL * (1 + 1j * ANOMALY_SIGMA_RATIO)

sigma = Scalar(sigma_bg, Function(dom))
sigma.setTaggedValue("Anomaly1", sigma_anomaly)
sigma.setTaggedValue("Anomaly2", sigma_anomaly)
sigma.setTaggedValue("Anomaly3", sigma_anomaly)

if not args.silent:
    print("sigma_bg = ", sigma_bg)
    print("sigma_anomaly = ", sigma_anomaly)
# ... SETUP SOLVER
solver = SIPSolver(dom, surface_mask=surface_mask,
                   rtol=1e-8, pdetol=1e-10, iter_max=100, verbose=not args.silent)
start_time = time.perf_counter()
# ... set the source potential
solver.setSource(source_location=SRC_LOC, tag="S1")
# .... set the primary potential
solver.setPrimaryPotential(sigma_bg, current=SRC_I)
# .... set the conductivity .....
solver.setConductivity(sigma, sigma_bc=sigma_bg)
end_time = time.perf_counter()
elapsed_time_setup = end_time - start_time
print(f"SIP Set-up time: {elapsed_time_setup:.4f} seconds (elapsed time)")

# .... solve it ...
start_time = time.perf_counter()
u = solver.solve()
end_time = time.perf_counter()
elapsed_time_solver = end_time - start_time
if not args.silent:
    print("u real : ", u.real())
    print("  imag : ", u.imag())
    print("history :", solver.history)
    print(f"SIP Solver time: {elapsed_time_solver:.4f} seconds (elapsed time)")

if args.test:
        pde = setupERTPDE(dom, tolerance=1e-10, isComplex=True)
        optionsG = pde.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        S = Scalar(0., DiracDeltaFunctions(dom))
        S.setTaggedValue("S1", SRC_I)
        pde.setValue(A=sigma * kronecker(3), d=sigma_bg * solver.alpha, y_dirac=S)
        start_time = time.perf_counter()
        du_test = pde.getSolution() - u
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        # ...............
        if not args.silent:
            print("Difference to direct solver results:")
            print(" real  = ", du_test.real())
            print(" imag  = ", du_test.imag())
            print(" total = ", Lsup(abs(du_test)) )
            print(f"Direct Solver time: {elapsed_time:.4f} seconds (elapsed time)")
#saveSilo("setup", m = anomaly_mask, u_src=u_src)

data = {
"mesh": args.meshfile,
"contrast": args.contrast,
"ratio": args.ratio,
"history" : solver.history,
"iterations" : len(solver.history)-1,
"timing_setup" : elapsed_time_setup,
"timing_solver" : elapsed_time_solver,
"timing_total" : elapsed_time_solver + elapsed_time_setup
}
json_string = json.dumps(data)
print(json_string)


