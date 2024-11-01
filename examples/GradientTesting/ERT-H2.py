#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

from wx.lib.agw.pyprogress import Continue

sys.path.append(os.getcwd())
from fingal import ERTInversionH2
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np


CONFIG="config-ERT-H2"
TABFN="ERT-H2.log"

import logging
from datetime import datetime

logger=logging.getLogger('ERT-H2')
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

costf = ERTInversionH2(domain, data=survey,
                       sigma_0_ref=config.sigma0_ref, reg_tol=None,
                       w1=config.regularization_w1, maskOuterFaces=mask_face, dataRTolDC=config.data_rtol,
                       pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                       useLogMisfitDC=config.use_log_misfit_DC, logger=logger)


tabfile=open(TABFN, 'w')
#for w1, with_misfit in [(0.,True ), (1., False), (1.e4, True)]:
for w1, with_misfit in [(1., False), (0., True), (1.e8, True)]:
    print("w1 = ", w1)
    costf.setW1(w1)
    costf.ignoreERTMisfit(not with_misfit)

    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    ppz=(x - inf(x))*(x - sup(x)) * (y - inf(y))* (y - sup(y))
    ppz/=sup(abs(ppz))
    ppx=(y - inf(y))*(y - sup(y)) *(z - inf(z))
    ppx/=sup(abs(ppx))
    ppy=(x - inf(x))*(x - sup(x)) *(z - inf(z))
    ppy/=sup(abs(ppy))
    #====
    r = length(domain.getX())
    M=Vector(0., ContinuousFunction(domain))
    M[0]= -0.2e-2* r/Lsup(r) * ppy *ppz
    M[1] = 0.1e-2 * r / Lsup(r) * ppx * ppz
    M[2] = 0.2e-2 * r / Lsup(r) * ppy * ppx

    for d in [ [ppy*ppz,0,0] , [0,ppx*ppz,0], [0,0,ppy*ppx], [ppy*ppz,ppx*ppz,ppy*ppx]  ] :
        dm=(x+y+0.5*z)/Lsup(r)/35
        dM = Vector(0., ContinuousFunction(domain))
        f="["
        for i, c in enumerate(d):
            dM[i] = dm * c
            if c is 0:
                f+="0 "
            else:
                f+="1 "
        f+="]"
        args=costf.getArgumentsAndCount(M)
        G=costf.getGradientAndCount(M, *args)
        Dex = costf.getDualProductAndCount(dM, G)
        J0=costf.getValueAndCount(M,  *args)
        print("J(m)=%e"%J0)
        print("gradient = %s"%str(G))
        b=[]
        t=[]
        tabfile.write(f".. w1 , = {w1}, with_misfit= {with_misfit}, d={f}  .............\n")
        tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)\n")
        for k in range(4, 13):
            a=0.5**k
            J=costf.getValueAndCount(M+a*dM)
            D=(J-J0)/a
            b.append(log(abs(D-Dex)))
            t.append(log(a))
            tabfile.write("%d      %e %e %e %e %e %e\n"%(k,J0, J, Dex, D, D-Dex, (D-Dex)/a) )
        m, c = np.linalg.lstsq(np.vstack([np.array(t), np.ones(len(t))]).T, b, rcond=None)[0]
        if m < 0.999:
            tabfile.write(f"WARNING: Poor convergence rate = {m}.\n")
        else:
            tabfile.write(f"Convergence rate = {m}.\n")
logger.info("All done - Have a nice day!")