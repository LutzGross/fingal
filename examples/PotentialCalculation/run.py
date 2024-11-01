#!/usr/bin/python3
from esys.escript import *
from esys.finley import ReadGmsh
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.weipa import saveVTK
from esys.escript.pdetools import Locator
import numpy as np

# line of electrodes to inspect for voltage:
LineOffset = 6
NumElectrodes =  64
ElectrodeSpacing = 3
LineLength = (NumElectrodes-1) * ElectrodeSpacing
# location of injection
X_A = [0.,0.,0]
I_A =1
# read mesh from gmsh:
domain=ReadGmsh('brick.msh', 3, diracPoints=[X_A], diracTags=['A'], optimize=True )

# overall conductivity
sigma0 = 0.0002
sigma = Scalar(sigma0, Function(domain))
# and higher condictivity in the anomaly:
sigma.setTaggedValue("anomaly", 0.02)

# set up PDE:
pde = LinearSinglePDE(domain, isComplex=False)
pde.setSymmetryOn()
optionsG = pde.getSolverOptions()
optionsG.setSolverMethod(SolverOptions.PCG)
optionsG.setTolerance(1e-8)
if hasFeature('trilinos'):
    optionsG.setPackage(SolverOptions.TRILINOS)
    optionsG.setPreconditioner(SolverOptions.AMG)
    #optionsG.setTrilinosParameter("problem:type", "Poisson-3D")
    optionsG.setTrilinosParameter("verbosity", "none")
    optionsG.setTrilinosParameter("number of equations", 1)
    optionsG.setTrilinosParameter("problem: symmetric", True)

# gradient of the primary potential:
X=Function(domain).getX()
gradu1=I_A/(4*np.pi*sigma0) * (X-X_A)/length(X-X_A)**3

# we fix the secondory potential on all faces with exception of the top:
bbx=getBoundingBox(domain)
x=domain.getX()
fix_u2 = whereZero(x[0]-bbx[0].min)+whereZero(x[0]-bbx[0].max)*\
  whereZero(x[1]-bbx[1].min)+whereZero(x[1]-bbx[1].max)+whereZero(x[2]-bbx[2].min)

# set the PDE coefficients:
pde.setValue(A = sigma * kronecker(3), X=(sigma0-sigma) * gradu1, q=fix_u2)
# and get the secondary potential:
u2=pde.getSolution()
print("secondary potential :", str(u2))

# write to VTK for 3D visualization:
saveVTK("results", u2=u2)

# and we want to plot the potential along the line of electrodes:
import matplotlib.pyplot as plt
grab_values=Locator(Solution(domain), [ (-LineLength/2 + k * ElectrodeSpacing ,   LineOffset, 0) for k in range(NumElectrodes)])

plt.clf()
plt.plot([X[0] for X in grab_values.getX()], grab_values(u2)  )
plt.title("Voltage due to conductivity anomaly")
plt.xlabel("offset [m]")
plt.ylabel("voltage [V]")
plt.savefig("profile.png")