#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

from leather import Scale

sys.path.append(os.getcwd())
from fingal import IP2InversionH1
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np
from esys.escript.pdetools import MaskFromBoundaryTag

CONFIG="config-IP2-H1_0"
TABFN="IP2-H1_0.log"
ESTTOL = 5.e-8
import logging
from datetime import datetime

logger=logging.getLogger('IP-H1_0')
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


costf = IP2InversionH1(domain, data=survey, sigma_0=Scalar(config.sigma0_ref, ContinuousFunction(domain)), Mn_ref=config.Mn_ref,
                       w1=config.regularization_w1IP,
                       maskZeroPotential=mask_face,
                       stationsFMT=config.stationsFMT, pde_tol=config.pde_tol,
                       useLogMisfitIP=config.use_log_misfit_IP,
                       dataRTolIP=config.data_rtol,
                       logclip=config.clip_property_function,
                       zero_mean_m=True,
                       logger=logger)

tabfile=open(TABFN, 'w')
# test the regularization first:
for w1, with_IPmisfit in [ (1, True), (0, False), (1e-2, False) ]:
    costf.setW1(w1)
    costf.ignoreIPMisfit(with_IPmisfit)
    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    pp=(x - inf(x))*(x - sup(x)) * (y - inf(y))* (y - sup(y))*(z - inf(z))
    #pp = (z - inf(z))
    pp/=sup(abs(pp))
    #====
    r= length(domain.getX())
    M=pp
    #print(abs(inner(grad(M[0]), grad(M[1])) / (length(grad(M[0]))*length(grad(M[1])))))
    #1/0
    ddm=(x+y+0.5*z)/Lsup(r)/3*10*pp/10000
    ddm = (x + y + 0.5 * z) / Lsup(r) * pp / 6
    args = costf.getArgumentsAndCount(M)
    G = costf.getGradientAndCount(M, *args)
    tabfile.write(
        f".. w1 , = {w1}, with_IPmisfit= {with_IPmisfit}  .............\n")
    dM=ddm
    Dex = costf.getDualProductAndCount(dM, G)
    J0=costf.getValueAndCount(M,  *args)
    print("J(m)=%e"%J0)
    print("gradient = %s"%str(G))
    b=[]
    x=[]

    tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)     rel. err. eval.\n")
    for k in range(4, 13):
        a=0.5**k
        J=costf.getValueAndCount(M+a*dM)
        D=(J-J0)/a
        if abs(J0-J) > ESTTOL * max(J0,J):
            b.append(log(abs(D-Dex)))
            x.append(log(a))
        tabfile.write("%d      %e %e %e %e %e %e %e\n"%(k,J0, J,Dex,  D, D-Dex, (D-Dex)/a, abs(J0-J)/max(J0,J)) )
    if len(x) > 0  :
        m, c = np.linalg.lstsq(np.vstack([np.array(x), np.ones(len(x))]).T, b, rcond=-1)[0]
        if m < 0.99:
            tabfile.write(f"WARNING: Poor convergence rate = {m} (# of data = {len(x)}).\n")
        else:
            tabfile.write(f"Convergence rate = {m} (# of data = {len(x)}).\n")
logger.info("All done - Have a nice day!")