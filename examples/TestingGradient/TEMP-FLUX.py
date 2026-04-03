#!/usr/bin/python3
from esys.escript import *
import importlib, os, sys

sys.path.append(os.getcwd())
from fingal import InversionIPByFluxWithWeakBC, ConductivityModelByTemperature
from fingal import readElectrodeLocations, readSurveyData, makeMaskForOuterSurface
from esys.finley import ReadMesh
import numpy as np
from esys.escript.pdetools import MaskFromBoundaryTag

CONFIG = "config-TEMP-FLUX"
TABFN = "TEMP-FLUX.log"
ESTTOL = 5e-8

import logging
from datetime import datetime

logger = logging.getLogger('TEMP-FLUX')
logger.setLevel(logging.DEBUG)
config = importlib.import_module(CONFIG)

elocations = readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
logger.info("%s electrode locations read from %s." % (len(elocations), config.stationfile))

domain = ReadMesh(config.meshfile)
logger.info("Mesh read from " + config.meshfile)

survey = readSurveyData(config.datafile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates,
                        columns=config.datacolumns,
                        dipoleInjections=config.dipoleInjections, dipoleMeasurements=config.dipoleMeasurements,
                        delimiter=config.datadelimiter, commend='#', printInfo=True)
assert survey.getNumObservations() > 0, "no data found."

mask_face = MaskFromBoundaryTag(domain, *config.faces_tags)

z = domain.getX()[2]
maskTopSurfaceElements = Scalar(0., FunctionOnBoundary(domain))
maskTopSurfaceElements.setTaggedValue('faces', 1)

costf = InversionIPByFluxWithWeakBC(domain, data=survey, maskZeroPotential=mask_face,
                                    conductivity_model=ConductivityModelByTemperature(),
                                    sigma_src=None,
                                    pde_tol=config.pde_tol, stationsFMT=config.stationsFMT,
                                    length_scale=config.regularization_length_scale,
                                    useLogMisfitDC=config.use_log_misfit_DC, dataRTolDC=config.data_rtol,
                                    useLogMisfitIP=config.use_log_misfit_IP, dataRTolIP=config.data_rtol,
                                    weightingMisfitDC=0.5,
                                    w1=config.regularization_w1DC, reg_tol=None,
                                    surface_flux=maskTopSurfaceElements*10,
                                    set_temperature=config.temperature_longrange,
                                    mask_set_temperature=whereZero(z-inf(z)),
                                    flux_variation_factor=config.surface_values_variation_factor,
                                    thermal_conductivity=config.thermal_conductivity,
                                    heat_production=0,
                                    logger=logger)

tabfile=open(TABFN, 'w')
for w1, with_ERTmisfit, with_IPmisfit in [(1., False, False), (0., True, False),  (0., False, True), (1.e-4, True, True)]:
    print("w1 = ", w1)
    costf.setW1(w1)
    costf.ignoreERTMisfit(not with_ERTmisfit)
    costf.ignoreIPMisfit(not with_IPmisfit )

    x = domain.getX()[0]
    y = domain.getX()[1]
    z = domain.getX()[2]
    # make sure that boundary conditions for m are met:
    ppz=(z - sup(z))  *(z - inf(z))
    ppz=(z - inf(z))
    ppz/=Lsup(abs(ppz))
    ppy=(y - inf(y))*(y - sup(y))
    ppy=(y - inf(y))
    ppy/=Lsup(abs(ppy))
    ppx=(x - inf(x))*(x - sup(x))
    ppx=(x - sup(x))

    ppx/=Lsup(abs(ppx))
    blob=abs(ppz*ppx*ppy)**8
    #====

    M=Vector(0., ContinuousFunction(domain))
    M[0]=  0.1
    M[1] = 0.2
    M[2] = -0.5 #* ppz

    #for d in [ [ppy*ppx,0,0] , [0,ppx*ppy,0], [0,0,ppz], [ppy*ppx,ppx*ppy,ppz]  ] :
    for d in [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]:
        dm = - blob*100
        dM = Vector(0., ContinuousFunction(domain))
        f="["
        for i, c in enumerate(d):
            if i < 2:
                fac=1
            else:
                fac =-1
            dM[i] = fac * dm * c
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
        tabfile.write(f".. w1 , = {w1}, with_ERTmisfit= {with_ERTmisfit}, with_IPmisfit ={with_IPmisfit}, d={f}  .............\n")
        tabfile.write("log(a)     J(m)        J(m+a*p)       grad        num. grad     error O(a)   O(1)\n")
        for k in range(4, 13):
            if k == 100:
                1/0
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
