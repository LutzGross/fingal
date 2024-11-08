#!/usr/bin/python3
from minetools import *
from esys.finley import ReadMesh
from esys.weipa import saveSilo
from esys.escript import *
from esys.escript.pdetools import Locator
from fingal import readElectrodeLocations, readSurveyData
import numpy as np
import os

import config
GEOFILE = "mine.geo"
PLOTDIR="plots"
SILODIR="results"
REL_ERR = 0.02

# extract some geometrical information from the geo file:
minegeo=getGeometryFromGeoFile(GEOFILE)
print(f"geometrical data read from {GEOFILE}")

elocations=readElectrodeLocations(config.stationfile, delimiter=config.stationdelimiter)
print("%s electrode locations read from %s."%(len(elocations), config.stationfile))

schedule = readSurveyData(config.schedulefile, stations=elocations, usesStationCoordinates=config.usesStationCoordinates,
                        columns=[],
                        dipoleInjections=config.dipoleInjections,
                        dipoleMeasurements=config.dipoleMeasurements,
                        delimiter=config.datadelimiter,
                        commend='#', printInfo=True)


domain=ReadMesh(config.meshfile)
print("Mesh read from "+config.meshfile)

grab_values_stations = Locator(DiracDeltaFunctions(domain),  [schedule.getStationLocationByKey(S) for S in schedule.getStationNumeration()])
# setup conductivity
RHO_BASE = 5000
RHO_COAL = 500
RHO_MINE = 5000

RHO_BASE = 500
RHO_COAL = 500
RHO_MINE = 500
RHO_GOAF = 50000
RHO_REF = sqrt(RHO_BASE*RHO_MINE)

rho = Scalar(RHO_REF, Function(domain) )
rho.setTaggedValue( 'Base', RHO_BASE)
rho.setTaggedValue('Seam', RHO_COAL)
rho.setTaggedValue( 'Goaf', RHO_GOAF)
rho.setTaggedValue('Mass', RHO_MINE)
rho.expand()
rho_raise_factor_damage = 5
rho, damage, damage_at_stations = applyDamage2(rho, minegeo, rho_raise_factor_damage = rho_raise_factor_damage, grab_values_stations=grab_values_stations)
sigma_at_stations = 1/(RHO_COAL *  (1 + damage_at_stations * (rho_raise_factor_damage-1) ) )


primary_potentials = makePrimaryPotentials(domain, minegeo, schedule=schedule, sigma_at_stations=sigma_at_stations)

secondary_potentials = makeSecondaryPotentials(domain, minegeo, sigma = 1/rho,
                                               primary_potentials=primary_potentials, schedule=schedule,  sigma_at_stations=sigma_at_stations)

potentials_at_stations = {}
for iA in primary_potentials:
    u1 = np.array(grab_values_stations( primary_potentials[iA]))
    u2 = np.array(grab_values_stations(secondary_potentials[iA]))
    potentials_at_stations[iA] = u1 + u2

f=open(config.datafile,'w')
for A,B,M,N in schedule.tokenIterator():
    iA=schedule.getStationNumber(A)
    iB=schedule.getStationNumber(B)
    iM=schedule.getStationNumber(M)
    iN=schedule.getStationNumber(N)
    if REL_ERR > 0 :
        pert = np.random.uniform(low=-REL_ERR, high=REL_ERR)
    else:
        pert = 0.
    u=(potentials_at_stations[iA][iM]-potentials_at_stations[iB][iM]-potentials_at_stations[iA][iN]+potentials_at_stations[iB][iN]) * (1+pert)
    f.write(f"{A:d}, {B:d}, {M:d}, {N:d}, {u:g}\n")
f.close()
u101=primary_potentials[schedule.getStationNumber(101)]+secondary_potentials[schedule.getStationNumber(101)]
u146=primary_potentials[schedule.getStationNumber(146)]+secondary_potentials[schedule.getStationNumber(146)]
u123=primary_potentials[schedule.getStationNumber(123)]+secondary_potentials[schedule.getStationNumber(123)]
saveSilo("setup", sigma=1/rho, u=u101-u146, u101=u101, u146=u146, u123=u123, damage=damage, tag=makeTagMap(Function(domain)))


