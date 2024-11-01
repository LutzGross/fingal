#!/usr/bin/python3
from esys.escript import *
from esys.escript.linearPDEs import LinearPDESystem, SolverOptions
from esys.finley import ReadMesh
from esys.weipa import saveSilo
import esys.escript.unitsSI as U
MESH="mesh.fly"
domain=ReadMesh(MESH)
print("Mesh read from "+MESH)


LAM=3.7*U.Giga * U.Pa
MU=8.7*U.Giga * U.Pa
Pmax=1 *U.Mega * U.Pa


mypde=LinearPDESystem(domain)
mypde.setSymmetryOn()
mypde.getSolverOptions().setVerbosityOn()
optionsG = mypde.getSolverOptions()
if hasFeature('trilinos'):
   optionsG.setPackage(SolverOptions.TRILINOS)
   #optionsG.setPreconditioner(SolverOptions.AMG)
   optionsG.setTrilinosParameter("verbosity", "none")
   optionsG.setTrilinosParameter("number of equations", 3)
   optionsG.setTrilinosParameter("problem: symmetric", True)

lam = Scalar(LAM, Function(domain))
lam.setTaggedValue('Seam', LAM/10)
lam.setTaggedValue('Mass', LAM)
lam.setTaggedValue( 'Goaf', LAM/10)
lam.setTaggedValue( 'Base', LAM)
mu = Scalar(MU, Function(domain))
mu.setTaggedValue('Seam', MU/10)
mu.setTaggedValue( 'Goaf', MU/10)
mu.setTaggedValue( 'Base', MU)
mu.setTaggedValue( 'Mass', MU)
#... set coefficients ...
C=mypde.createCoefficient('A')
for i in range(domain.getDim()):
  for j in range(domain.getDim()):
     C[i,i,j,j]+=lam
     C[i,j,i,j]+=mu
     C[i,j,j,i]+=mu
rho = Scalar(2500, Function(domain))

x=domain.getX()
#msk = whereZero(x[2]-inf(x[2]))*[1.,1.,1.]
msk= whereZero(x[0]-inf(x[0])) *[1.,0.,0.] \
   + whereZero(x[1]-inf(x[1])) *[0.,1.,0.] \
   +  whereZero(x[2]-inf(x[2]))*[0.,0.,1.]

n=domain.getNormal()
X=n.getFunctionSpace().getX()
P=-Pmax * wherePositive(whereZero(X[0]-sup(x[0]))+0*whereZero(X[1]-sup(x[1])))- Pmax*0 * whereZero(X[2]-sup(x[2])) \
+Pmax*wherePositive(whereZero(X[0]-inf(x[0]))+0*whereZero(X[1]-inf(x[1])))
mypde.setValue(A=C,y=P * n,q=msk, Y=-rho*9.81*0* [0,0,1])
#... solve pde ...
u=mypde.getSolution()
# stress:
g=grad(u, where=ReducedFunction(domain))
stress = 2 * mu * symmetric(g) + lam * trace(g) * kronecker(3)

q=0.3
K=30000
S = sqrt(1/2)*length(deviatoric(stress)) - q* trace(stress)/3
#S= sqrt( 0.5*(stress[0,0]-stress[1,1])**2 + 0.5*(stress[0,0]-stress[2,2])**2 + 0.5*(stress[2,2]-stress[1,1])**2 + 3*stress[1,0]**2 + 3*stress[1,2]**2 + 3*stress[2,0]**2  )
#... output ...
saveSilo("deform.silo",u=u, vanMis=S)