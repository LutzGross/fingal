#!/usr/bin/python3

from esys.finley import ReadMesh
from esys.weipa import saveSilo
from esys.escript import *
from esys.escript.pdetools import Locator
from fingal import readElectrodeLocations, readSurveyData, setupERTPDE
import numpy as np
import os

import config
RAISE_FACTOR_DAMAGED_RESISTIVITY = 5
#GEOFILE = "mine.geo"
PLOTDIR="plots"
SILODIR="results"
REL_ERR = 0.02

# extract some geometrical information from the geo file:
#minegeo=getGeometryFromGeoFile(GEOFILE)
#print(f"geometrical data read from {GEOFILE}")

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

# setup base resistivity:
RHO_REF = 100
RHO_BASE = RHO_REF
RHO_COAL = RHO_REF * 2
RHO_MINE = RHO_REF
RHO_GOAF = RHO_REF * 50


rho = Scalar(RHO_REF, Function(domain) )
rho.setTaggedValue( 'Base', RHO_BASE)
rho.setTaggedValue('Seam', RHO_COAL)
rho.setTaggedValue( 'Goaf', RHO_GOAF)
rho.setTaggedValue('Core', RHO_MINE)
rho.expand()



## .... set damage ............................
## all in meter
# offset of the Damaged zone from the SeamEnd
DamageZoneOffset  = 100
DamageHeight = 80
DamageGeometryExponent = 0.5
DamageBaseDepth = 40;
FringeWidth = 10;
ResetDamagedZoneSouth = 3 * FringeWidth
ResetDamagedZoneNorth = 3* FringeWidth
# ========
SeamHeight = 13
SeamStart = -115.23432427900843 # Point 1199
SeamEnd = 147.48472712998046 # Point 314
SeamNorth = 56.82856875646394 # # Point 1214
SeamSouth = -143.17124071449507 # Point 1199
MinedEnd = 197.48471409198828 # Point 404
GoafEnd = 241.06384326395346 # Point 462


#Point(1199) = {-115.23432427900843, -143.17110019456595, SeamHeight, meshSizeSeam};
#Point(314) = {147.48472712998046, -143.17124071449507, 0., meshSizeSeam};
#Point(404) = {197.48471409198828, -143.17125238850713, 0., meshSizeCore};
#Point(462) = {241.06384326395346, -143.17125472647604, 0., meshSizeCore};
#Point(1214) = {-115.23432427900843, 56.82856875646394, 0., meshSizeSeam};

X = ReducedFunction(domain).getX()
x = X[0]
y = X[1]
z = X[2]

#GoafEnd = minegeo.UnMinedLength + minegeo.ExtractionLength + minegeo.GoafLength
DamagedZoneStart = SeamEnd - DamageZoneOffset

h_top = DamageHeight * clip((x - DamagedZoneStart) / (GoafEnd - DamagedZoneStart), minval=0.) ** DamageGeometryExponent
h_base = DamageBaseDepth * clip((x - DamagedZoneStart) / (GoafEnd - DamagedZoneStart), minval=0.) ** DamageGeometryExponent
m2 = whereNonNegative(X[2] + h_base) * whereNonPositive(X[2] - h_top)
m2 *= whereNonPositive(X[1] - (SeamNorth-ResetDamagedZoneNorth) ) * whereNonNegative(X[1] -(SeamSouth+ResetDamagedZoneSouth))

# no damage in the goaf:
#m2.setTaggedValue('Goaf', 0)
#m2.setTaggedValue('Padding', 0)

# apply some smoothing:
pde = setupERTPDE(domain, tolerance=1e-8)
pde.setValue(D=1, A=FringeWidth ** 2 * kronecker(3), Y=interpolate(m2, Function(domain)))
m = clip(pde.getSolution(), minval=0, maxval=1)

m3 = clip(interpolate(m, rho.getFunctionSpace()), minval=0, maxval=1)
m3.setTaggedValue('Goaf', 0)
m3.setTaggedValue('Padding', 0)

rho = rho * (1 + m3 * (RAISE_FACTOR_DAMAGED_RESISTIVITY - 1))


rho_stations = Scalar(RHO_COAL, grab_values_stations.getFunctionSpace())
rho_stations +=  whereNegative(grab_values_stations.getFunctionSpace().getX()[1]-(SeamSouth-1))  * (RHO_MINE-rho_stations)
rho_stations = rho_stations * (1 + interpolate(m, grab_values_stations.getFunctionSpace()) * (RAISE_FACTOR_DAMAGED_RESISTIVITY - 1))

rho_src=RHO_COAL

#m3, np.array(grab_values_stations(m))
saveSilo("setup", sigma=1/rho, tag=makeTagMap(Function(domain)))

#saveSilo("X", rho=rho, m3=m)
1/0


src_potentials = getSourcePotentials(domain,
                                     sigma=1/rho_src,
                                     survey=schedule,
                                     maskZeroPotential=q,
                                     stationsFMT='s%d')
pde = setupERTPDE(domain, tolerance=1e-10)
pde.setValue(A=sigma * kronecker(3), q=q)
potential, src_potential_scale = getAdditivePotentials(pde, sigma=1/rho,
                                                        schedule=schedule,
                                                        sigma_stations=1/rho_stations,
                                                        src_potentials=src_potentials,
                                                        sigma_src=1/rho_src)
#==============================
1/0
rho, damage, damage_at_stations = applyDamage2(rho, minegeo, rho_raise_factor_damage = rho_raise_factor_damage, grab_values_stations=grab_values_stations)
sigma_at_stations = 1/(RHO_COAL *  (1 + damage_at_stations * (rho_raise_factor_damage-1) ) )


primary_potentials = makePrimaryPotentials(domain, minegeo, schedule=schedule, sigma_at_stations=sigma_at_stations)

secondary_potentials = makeSecondaryPotentials(domain, minegeo, sigma = 1/rho, sigma_faces = 1/RHO_REF,\
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


