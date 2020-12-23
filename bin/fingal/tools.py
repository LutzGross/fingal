"""
Fingal - some tools 

by l.gross@uq.edu.au, Dec 2020.
"""


from esys.escript import Scalar

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
