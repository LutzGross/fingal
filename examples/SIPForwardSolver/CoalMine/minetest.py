#!/usr/bin/python3
from esys.escript import *
from esys.finley import ReadMesh
from esys.weipa import saveSilo
from fingal import setupERTPDE, SIPSolver
from esys.escript.linearPDEs import SolverOptions
import time, json


MESHFIILE = "./mine.fly"
SRC_LOC = (44.1925003270153, 56.8285531244474, 1.5)

dom = ReadMesh(MESHFIILE)
print("mesh read from ", MESHFIILE)

# ... mask of the insulated surface.
surface_mask = Scalar(0., FunctionOnBoundary(dom))
# ... set-up conductivity distribution:
# ...... the overall conductiviry also acts as primary conductivity
sigma_bg = 1e-3 * (2.5 + 1j * 0.)
# ....... conductivity in the anomaly

sigma_coal = 1e-3 * (0.5 + 1j * 0.15)

sigma = Scalar(sigma_bg, Function(dom))
sigma.setTaggedValue("Seam", sigma_coal)

print("sigma_bg = ", sigma_bg)
print("sigma_coal = ", sigma_coal)




# ... SETUP SOLVER
solver = SIPSolver(dom, surface_mask=surface_mask,
                   rtol=1e-8, pdetol=1e-10, iter_max=100, verbose=True)
start_time = time.perf_counter()
# ... set the source potential
solver.setSource(source_location=SRC_LOC, tag="S1")
# .... set the primary potential
solver.setPrimaryPotential(sigma_bg, current=1)
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
print("u real : ", u.real())
print("  imag : ", u.imag())
print("history :", solver.history)
print(f"SIP Solver time: {elapsed_time_solver:.4f} seconds (elapsed time)")

if False:
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
        print("Difference to direct solver results:")
        print(" real  = ", du_test.real())
        print(" imag  = ", du_test.imag())
        print(" total = ", Lsup(abs(du_test)) )
        print("Direct Solver time: {elapsed_time:.4f} seconds (elapsed time)")
#saveSilo("setup", m = anomaly_mask, u_src=u_src)

data = {
"mesh": MESHFIILE,
"history" : solver.history,
"iterations" : len(solver.history)-1,
"timing_setup" : elapsed_time_setup,
"timing_solver" : elapsed_time_solver,
"timing_total" : elapsed_time_solver + elapsed_time_setup
}
json_string = json.dumps(data)
print(json_string)


