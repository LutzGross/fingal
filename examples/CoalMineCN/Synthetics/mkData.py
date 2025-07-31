#!/usr/bin/python3

from esys.finley import ReadMesh
from esys.weipa import saveSilo
from esys.escript import *
from esys.escript.pdetools import Locator
from fingal import readElectrodeLocations, readSurveyData, setupERTPDE, getSourcePotentials, getAdditivePotentials
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
RHO_REF = 200
RHO_BASE = RHO_REF
RHO_COAL = RHO_REF
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
DamageZoneOffset  = 140 #  220

DamageHeight = 60
DamageGeometryExponent = 0.5
DamageBaseDepth = 20;

FringeWidth = 5;
ResetDamagedZoneSouth = 2 * FringeWidth
ResetDamagedZoneNorth = 2* FringeWidth
FringeWidth = 1;
ResetDamagedZoneSouth = 0 * FringeWidth
ResetDamagedZoneNorth = 0* FringeWidth
#ResetDamagedZoneSouth = 100
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


rho_src_stations = Scalar(RHO_COAL, grab_values_stations.getFunctionSpace())
rho_src_stations +=  whereNegative(grab_values_stations.getFunctionSpace().getX()[1]-(SeamSouth-1))  * (RHO_MINE-rho_src_stations)
rho_stations = rho_src_stations * (1 + interpolate(m, grab_values_stations.getFunctionSpace()) * (RAISE_FACTOR_DAMAGED_RESISTIVITY - 1))
rho_src=RHO_COAL

#m3, np.array(grab_values_stations(m))


#saveSilo("X", rho=rho, m3=m)
x=domain.getX()
q=Scalar(0., x.getFunctionSpace())
for i in range(3):
    q+=whereZero(x[i]-inf(x[i]))+whereZero(x[i]-sup(x[i]))

src_potentials = getSourcePotentials(domain,
                                     sigma=1/rho_src,
                                     survey=schedule,
                                     maskZeroPotential=q,
                                     stationsFMT='s%d')
print(f"{len(src_potentials)} source potentials calculated.")

pde = setupERTPDE(domain, tolerance=1e-10)
pde.setValue(A=1/rho * kronecker(3), q=q)
potential2, src_potential_scale = getAdditivePotentials(pde, sigma=1/rho,
                                                        schedule=schedule,
                                                        sigma_stations=grab_values_stations(1/rho_stations),
                                                        source_potential=src_potentials,
                                                        sigma_src=1/rho_src,
                                                        sigma_src_stations=grab_values_stations(1/rho_src_stations))

print(f"{len(potential2)} secondary potentials calculated.")


potentials_at_stations = {}
for iA in src_potentials:
    u1 = np.array(grab_values_stations( src_potentials[iA]))
    u2 = np.array(grab_values_stations(potential2[iA]))
    a = src_potential_scale[iA]
    potentials_at_stations[iA] = a * u1 + u2

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
print(f"data saved to file {config.datafile}.")

kargs = { 'sigma' : 1/rho,   'damage' : m3, 'tag' : makeTagMap(Function(domain)) }
for p in [] : # [101, 140, 220] :
    kargs[f'u{p}'] =   src_potentials[schedule.getStationNumber(p)]
#                    * src_potential_scale[schedule.getStationNumber(p)] + \
#                      potential2[schedule.getStationNumber(p)]
    print(f"potential u{p} added to output. ")
saveSilo("setup", **kargs)


