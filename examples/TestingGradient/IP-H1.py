#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

sys.path.append(os.getcwd())
from fingal import IPInversionH1
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np
from esys.escript.pdetools import MaskFromBoundaryTag

CONFIG="config-IP-H1"
TABFN="IP-H1.log"

import logging
from datetime import datetime

logger=logging.getLogger('IP-H1')
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

mask_face = MaskFromBoundaryTag(domain, *config.faces_tags)

costf=IPInversionH1(domain, data=survey,
                         sigma_0_ref=config.sigma0_ref, Mn_ref = config.Mn_ref,
                         w1=config.regularization_w1,
                        maskZeroPotential=mask_face, dataRTolDC= config.data_rtol,
                         pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                         useLogMisfitDC= config.use_log_misfit_DC, logger=logger)

tabfile=open(TABFN, 'w')
#for w1, with_misfit in [(0.,True ), (1., False), (1.e4, True)]:
for w1, theta, with_ERTmisfit, with_IPmisfit  in [(0., 0., True, False), (0., 0., False, True), (1., 0., False, False), (0., 1., False, False) ]:
#for w1, theta, with_ERTmisfit, with_IPmisfit in [(0., 0., False, True)]:
    print("w1 = ", w1)

    costf.setW1andTheta(w1, theta)
    costf.ignoreERTMisfit(not with_ERTmisfit)
    costf.ignoreIPMisfit(not with_IPmisfit)

    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    pp=(x - inf(x))*(x - sup(x)) * (y - inf(y))* (y - sup(y))*(z - inf(z))
    #pp = (z - inf(z))
    pp/=sup(abs(pp))
    #====
    r= length(domain.getX())
    M=RandomData((2,), ContinuousFunction(domain))*pp
    #print(abs(inner(grad(M[0]), grad(M[1])) / (length(grad(M[0]))*length(grad(M[1])))))
    #1/0
    ddm=(x+y+0.5*z)/Lsup(r)/3*10*pp/10000
    ddm = (x + y + 0.5 * z) / Lsup(r) * pp / 6
    for d in [ [1, 0] , [0, -0.25] , [1, -0.5] ] :
        tabfile.write(
            f".. w1 , = {w1}, theta = {theta}, with_ERTmisfit= {with_ERTmisfit}, with_IPmisfit= {with_IPmisfit} d={d}  .............\n")
        dM=ddm * d
        args=costf.getArgumentsAndCount(M)
        G=costf.getGradientAndCount(M, *args)
        Dex = costf.getDualProductAndCount(dM, G)
        J0=costf.getValueAndCount(M,  *args)
        print("J(m)=%e"%J0)
        print("gradient = %s"%str(G))
        b=[]
        x=[]
        please_fit = True
        tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)\n")
        for k in range(4, 13):
            a=0.5**k
            J=costf.getValueAndCount(M+a*dM)
            D=(J-J0)/a
            if abs(D-Dex) > 0:
                b.append(log(abs(D-Dex)))
                x.append(log(a))
            else:
                please_fit = False
            tabfile.write("%d      %e %e %e %e %e %e\n"%(k,J0, J,Dex,  D, D-Dex, (D-Dex)/a) )
        if please_fit :
            m, c = np.linalg.lstsq(np.vstack([np.array(x), np.ones(len(x))]).T, b, rcond=-1)[0]
            if m < 0.999:
                tabfile.write(f"WARNING: Poor convergence rate = {m}.\n")
            else:
                tabfile.write(f"Convergence rate = {m}.\n")
logger.info("All done - Have a nice day!")