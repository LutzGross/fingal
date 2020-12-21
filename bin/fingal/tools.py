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
