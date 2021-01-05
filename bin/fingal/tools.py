"""
Fingal - some tools 

by l.gross@uq.edu.au, Dec 2020.
"""


from esys.escript import Scalar, getMPIRankWorld, integrate, hasFeature
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions

def makeTagField(functionspace):
    """
    this creates a data object making the tagged regions with the corresponding tag id 
    
    :param functionspace: function space to create tags for
    :type functionspace: `esys.escript.FunctionSpace` 
    :return: `esys.escript.Scalar` with tags 
    """
    out=Scalar(-1,functionspace)
    for t in functionspace.getListOfTags():
        out.setTaggedValue(t,t)
    out.expand()
    return out

def FindNearestElectrode(x, y, z, electrodes={}):
    """
    finds the nearest electrode in the dictionary electrodes
    (x,y,z) is the coordinates of a proposed electrode
    :param x: x coordinate of electrode
    :param y: y coordinate of electrode
    :param z: x coordinate of electrode
    :param electrodes: dictionary of defined electrodes: electrodes[ide]=X=(xe, ye, ze)
    :return: id of nearest electrode and distance and distance to 
    
    """
    distmin=1e88
    idemin=None
    for ide in electrodes:
        X=electrodes[ide]
        dist=((X[0]-x)**2 + (X[1]-y)**2 + (X[2]-z)**2)**0.5
        if dist < distmin :
           distmin=dist
           idemin=ide
    if not idemin is None:
        return int(idemin), distmin
    else:
        return idemin, distmin


def setupERTPDE(domain, tolerance=1e-8, poisson=True, debug=0):
    """
    used to setup all PDEs fro ERT related inversion. If available TRILINOS is usered.
    
    :param domain: domain of the PDE
    :type domain: `esys.escript.AbstractDomain`
    :param tolerance: solver tolerance
    :type tolerance: `float` 
    :param poisson: if True TRILINOS settings fro Poisson problems is set.
    :type poisson: `bool`
    :return: linear, scalar PDE with real coefficients.
    :rtype: `esys.escript.linearPDEs.LinearPDE`
    

    """
    pde=LinearSinglePDE(domain, isComplex=False)
    pde.setSymmetryOn()
    optionsG=pde.getSolverOptions()
    optionsG.setSolverMethod(SolverOptions.PCG)
    optionsG.setTolerance(tolerance)
    if hasFeature('trilinos'):
        if debug and getMPIRankWorld() == 0:
            print("TRILINOS solver used.")
        optionsG.setPackage(SolverOptions.TRILINOS)
        optionsG.setPreconditioner(SolverOptions.AMG)
        if poisson:
            optionsG.setTrilinosParameter("problem:type", "Poisson-3D")
        optionsG.setTrilinosParameter("verbosity", "none")
        optionsG.setTrilinosParameter("number of equations", 1)
        optionsG.setTrilinosParameter("problem: symmetric", True)
 
    return pde


def getR2(y, y_fitted, chi=None):
        """
        calculates the coefficient of determination R^2 for  `y_fitted` as prediction for `y` over a region marked by chi>0 defined by
        
        R^2=1 - S_res/S_tot
        
        with S_res=int(chi*(y-y_fitted*1)**2,  S_tot=int(chi*(y-m(y)*1)**2), m(y)=int(chi*y)/int(chi)
        
        If R^2=1 then `y_fitted` is predicts `y` exactly. If R^2 then `y_fitted` does not make a better prediction than the mean. 
        
        :param y: target distribution
        :type y: `esys.escript.Scalar`
        :param y_fitted: fitted distribution
        :type y_fitted: `esys.escript.Scalar`
        :param chi: marker/weighting for region of interest
        :type chi: `esys.escript.Scalar` or None
        :rtype: `float`
        """
        if chi is None:
            chi=Scalar(1., Function(y_fitted.getFunctionSpace().getDomain()))
    
        ybar=integrate(chi*y)/integrate(chi)
        S_res=integrate(chi*(y-y_fitted)**2)
        S_tot=integrate(chi*(y-ybar)**2)
        if S_tot > 0:
            R2=1-S_res/S_tot
        else:
            if S_res > 0:
                R2=0.
            else:
                R2=1.
            
        return R2
