"""
this interpolates to a given domain
"""
from esys.escript import *
import numpy as np
from math import floor
from scipy.interpolate import RegularGridInterpolator

def mapToDomain(domain, data_table, resolution=1., origin=(0.,0., 0.), data0=0., locators=None, tags0=[]):
    """
    this maps the data_table values into the domain into the data object data.
    the table is assumed to be cell centered.
    """
    ndim=len(data_table.shape)
    points = tuple( [ np.linspace(origin[i], origin[i]+N*resolution, N) for i,N in enumerate(data_table.shape) ])
    interpolator=RegularGridInterpolator(points, data_table, method='linear', bounds_error=False, fill_value=data0)

    X=domain.getX()
    if ndim == 3:
        xi=[X.getTupleForDataPoint(p) for p in range(X.getNumberOfDataPoints()) ]
    else:
        xi=[X.getTupleForDataPoint(p)[::2] for p in range(X.getNumberOfDataPoints()) ]

    datai=interpolator(xi)
    data1=Scalar(data0, X.getFunctionSpace())
    data1.expand()
    [  data1.setValueOfDataPoint(p,datai[p]) for p in range(data1.getNumberOfDataPoints()) ]
    if locators:
        out=locators(data1)
    else:
        out=[]
    data=interpolate(data1, Function(domain))
    [ data.setTaggedValue(t , data0) for t in tags0 ]
    
    return data, out
