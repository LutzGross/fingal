#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

sys.path.append(os.getcwd())
from fingal import ERTInversionGauss, ERTInversionGaussWithDiagonalHessian
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
from esys.escript.pdetools import MaskFromBoundaryTag

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

mask_face = MaskFromBoundaryTag(domain, *config.faces_tags)

#costf=ERTInversionGaussWithDiagonalHessian(
costf=ERTInversionGauss(
                        domain, data=survey,
                         sigma_0_ref=config.sigma0_ref,
                         w1=config.regularization_w1DC, penalty_factor=config.regularization_penalty_factor,
                         maskZeroPotential= mask_face, dataRTolDC= config.data_rtol,
                         pde_tol=config.pde_tol, stationsFMT=config.stationsFMT, logclip=config.clip_property_function,
                         useLogMisfitDC= config.use_log_misfit_DC, logger=logger)

tabfile=open(TABFN, 'w')
for w1, with_misfit in [(0.,True ), (1., False), (1.e4, True)]:
    print("w1 = ", w1)
    costf.setW1(w1)
    costf.ignoreERTMisfit(not with_misfit)

    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    ppz = (x - inf(x)) * (x - sup(x)) * (y - inf(y)) * (y - sup(y))
    ppz /= sup(abs(ppz))
    ppx = (y - inf(y)) * (y - sup(y)) * (z - inf(z))
    ppx /= sup(abs(ppx))
    ppy = (x - inf(x)) * (x - sup(x)) * (z - inf(z))
    ppy /= sup(abs(ppy))
    # ====
    r = length(domain.getX())
    M = Data(0., (4,), ContinuousFunction(domain))
    M[0] = 0.67e-2 * r / Lsup(r) * ppy * ppz * ppx
    M[1] = -0.2e-2 * r / Lsup(r) * ppy * ppz
    M[2] = 0.1e-2 * r / Lsup(r) * ppx * ppz
    M[3] = 0.2e-2 * r / Lsup(r) * ppy * ppx

    for d in [[ppy * ppz * ppx, 0, 0, 0], [0, ppy * ppz, 0, 0], [0, 0, ppx * ppz, 0], [0, 0, 0, ppy * ppx],
              [ppy * ppz * ppx, ppy * ppz, ppx * ppz, ppy * ppx]]:
        dm = (x + y + 0.5 * z) / Lsup(r) / 40/10
        dM = Data(0., (4,), ContinuousFunction(domain))
        f = "["
        for i, c in enumerate(d):
            dM[i] = dm * c
            if c is 0:
                f += "0 "
            else:
                f += "1 "
        f += "]"
        args = costf.getArgumentsAndCount(M)
        G = costf.getGradientAndCount(M, *args)
        Dex = costf.getDualProductAndCount(dM, G)
        J0 = costf.getValueAndCount(M, *args)
        print("J(m)=%e" % J0)
        print("gradient = %s" % str(G))
        please_fit = True
        b = []
        t = []
        tabfile.write(f".. w1 , = {w1}, with_misfit= {with_misfit}, d={f}  .............\n")
        tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)\n")
        for k in range(4, 13):
            a = 0.5 ** k
            J = costf.getValueAndCount(M + a * dM)
            D = (J - J0) / a
            if abs(D - Dex) > 0:
                b.append(log(abs(D - Dex)))
                t.append(log(a))
            else:
                please_fit = False
            tabfile.write("%d      %e %e %e %e %e %e\n" % (k, J0, J, Dex, D, D - Dex, (D - Dex) / a))
        if please_fit:
            m, c = np.linalg.lstsq(np.vstack([np.array(t), np.ones(len(t))]).T, b, rcond=None)[0]
            if m < 0.999:
                tabfile.write(f"WARNING: Poor convergence rate = {m}.\n")
            else:
                tabfile.write(f"Convergence rate = {m}.\n")
logger.info("All done - Have a nice day!")