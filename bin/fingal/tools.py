from esys.escript import Scalar

def makeTagField(functionspace):
    """
    this creates a data object making the tagges regions with the corresponding tag id 
    """
    out=Scalar(-1,functionspace)
    for t in functionspace.getListOfTags():
        out.setTaggedValue(t,t)
    out.expand()
    return out

def FindNearestElectrode(x, y, z, electrodes):
    """
    (x,y,z) is the coordinates of a proposed electrode (by true data file)
    electrodes dictonary of defined electrodes in mesh electrodes[ide]=X=(xe,, ye, ze)
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
