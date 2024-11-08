#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

sys.path.append(os.getcwd())
from fingal import ERTInversionGauss
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np


CONFIG="config-ERT-Gauss"
TABFN="ERT-Gauss.log"

import logging
from datetime import datetime

logger=logging.getLogger('ERT-Gauss')
logger.setLevel(logging.DEBUG)
config = importlib.import_module(CONFIG)


elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
logger.info("%s electrode locations read from %s."%(len(elocations), config.stationfile))

domain=ReadMesh(config.meshfile)
logger.info("Mesh read from "+config.meshfile)


survey=readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates, columns=config.datacolumns,
                     dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements,
                      delimiter=config.datadelimiter, commend='#', printInfo=True)
assert survey.getNumObservations()>0, "no data found."

mask_face=makeMaskForOuterSurface(domain, taglist=config.faces_tags)

costf=ERTInversionGauss(domain, data=survey,
                         sigma_0_ref=config.sigma0_ref,
                         w1=config.regularization_w1,
                         maskOuterFaces= mask_face, dataRTolDC= config.data_rtol,
                         pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                         useLogMisfitDC= config.use_log_misfit_DC, logger=logger)

tabfile=open(TABFN, 'w')
for w1, with_misfit in [(0.,True ), (1., False), (1.e1, True)]:
    print("w1 = ", w1)
    tabfile.write(f".. w1 , = {w1}, with_misfit= {with_misfit} .............\n")
    costf.setW1(w1)
    costf.ignoreERTMisfit(not with_misfit)

    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    pp=(x - inf(x))*(x - sup(x)) * (y - inf(y))* (y - sup(y))*(z - inf(z))
    #pp = (z - inf(z))
    pp/=sup(abs(pp))
    #====

    r=length(domain.getX())
    m=r/Lsup(r)*pp

    ddm=(x+y+0.5*z)/Lsup(r)/3*pp
    dm=ddm

    args=costf.getArgumentsAndCount(m)
    G=costf.getGradientAndCount(m, *args)
    Dex = costf.getDualProductAndCount(dm, G)
    J0=costf.getValueAndCount(m,  *args)
    print("J(m)=%e"%J0)
    print("gradient = %s"%str(G))
    b=[]
    x=[]
    tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)\n")
    for k in range(4, 13):
        a=0.5**k
        J=costf.getValueAndCount(m+a*dm)
        D=(J-J0)/a
        b.append(log(abs(D-Dex)))
        x.append(log(a))
        tabfile.write("%d      %e %e %e %e %e %e\n"%(k,J0, J,Dex,  D, D-Dex, (D-Dex)/a) )
    m, c = np.linalg.lstsq(np.vstack([np.array(x), np.ones(len(x))]).T, b, rcond=1)[0]
    if m < 0.999:
        tabfile.write(f"WARNING: Poor convergence rate = {m}.\n")
    else:
        tabfile.write(f"Convergence rate = {m}.\n")
logger.info("All done - Have a nice day!")