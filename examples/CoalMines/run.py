#!/usr/bin/python3
from minetools import *
from esys.finley import ReadMesh
from esys.weipa import saveSilo
from esys.escript import *
from esys.escript.pdetools import Locator
import numpy as np
import os

GEOFILE = "mine.geo"
FLYFILE = "mine.fly"
PLOTDIR="plots"
SILODIR="results"

survey = { (101, 120) : 5 }
rho_ref = 100
rho_frac = rho_ref*100
# extract some geometrical information from the geo file:
minegeo=getGeometryFromGeoFile(GEOFILE)

minegeo.ResetFactureZoneNorth = 20
minegeo.ResetFactureZoneSouth = 20
minegeo.ThicknessFactureZone = 5
minegeo.FringFractureZone = 10
SteppingFractureZone = 30

# get the mesh
domain = ReadMesh(FLYFILE)
print(f"read FEM file from {FLYFILE}.")

north_grab_values_stations = Locator(DiracDeltaFunctions(domain), [ minegeo.Stations[k] for k in minegeo.LineNorth ])
south_grab_values_stations = Locator(DiracDeltaFunctions(domain), [ minegeo.Stations[k] for k in minegeo.LineSouth ])

# calculate primary potentials for constant conductivity sigma_ref
primary_potentials = makePrimaryPotentials(domain, minegeo, sigma_ref=1/rho_ref, survey=survey)
# ---------------- this is a bit of a hack: ----------------------------
A, B=list(survey.keys())[0]
I = survey[ (A,B) ]

Nsteps = 6
u1 = (primary_potentials[A] - primary_potentials[B]) * I
u1_north = north_grab_values_stations(u1)
u1_south = south_grab_values_stations(u1)
X_north, data1_north = makeMeasurements(minegeo.LineNorth, u1_north, minegeo, injections=[A, B], dir=0)
X_south, data1_south = makeMeasurements(minegeo.LineSouth, u1_south, minegeo, injections=[A, B], dir=0)
X_north, rho1_north= makeApparentResitivity(minegeo.LineNorth, u1_north, minegeo, injections=(A, B), I=I, dir=0)
X_south, rho1_south= makeApparentResitivity(minegeo.LineSouth, u1_south, minegeo, injections=(A, B), I=I, dir=0)
data2_north = {}
data2_south = {}
data_north = {}
data_south = {}
rho_north = {}
rho_south = {}
for kk in range(Nsteps):
    rho = makeResistivity1(domain, SteppingFractureZone * kk, rho_ref, rho_frac, minegeo)
    #===========
    secondary_potentials = makeSecondaryPotentials(domain, minegeo, sigma = 1/rho, sigma_ref=1/rho_ref, primary_potentials=primary_potentials)

    u2=(secondary_potentials[A]-secondary_potentials[B])*I
    u=u1+u2

    u2_north=north_grab_values_stations(u2)
    u2_south=south_grab_values_stations(u2)
    u_north=north_grab_values_stations(u)
    u_south=south_grab_values_stations(u)

    X_north, data2_north[kk] = makeMeasurements(minegeo.LineNorth , u2_north, minegeo, injections=[A, B], dir=0)
    X_south, data2_south[kk] = makeMeasurements(minegeo.LineSouth , u2_south, minegeo, injections=[A, B], dir=0)
    X_north, data_north[kk] = makeMeasurements(minegeo.LineNorth , u_north, minegeo, injections=[A, B], dir=0)
    X_south, data_south[kk] = makeMeasurements(minegeo.LineSouth , u_south, minegeo, injections=[A, B], dir=0)
    X_north, rho_north[kk] = makeApparentResitivity(minegeo.LineNorth, u_north, minegeo, injections=(A, B), I = I, dir=0)
    X_south, rho_south[kk] = makeApparentResitivity(minegeo.LineSouth,u_south, minegeo, injections=(A, B), I = I, dir=0)

    SILOFILE=os.path.join(SILODIR,f"u{kk}")
    saveSilo(SILOFILE, u1=u1, u=u, u2=u2,  rho=rho, tag=makeTagMap(ReducedFunction(domain)))
    print(f'results written to file {SILOFILE}.silo.')
import matplotlib.pyplot as plt


plt.scatter(X_north, np.array(rho1_north), label = "i" , s= 10)
for kk in rho_north:
    plt.scatter(X_north, np.array(rho_north[kk]), label = f"f{kk}", s= 10 )
plt.title(f"apparent resistivity Northern Line A,B={A},{B}")
plt.legend()
plt.xlabel("offset [m]")
plt.ylabel("rho [Ohmm]")
plt.savefig(os.path.join(PLOTDIR,"rho_north.png"))

plt.clf()
plt.scatter(X_south, np.array(rho1_south), label = "i", s= 10 )
for kk in rho_south:
    plt.scatter(X_south, np.array(rho_south[kk]), label = f"f{kk}", s= 10 )
plt.title(f"apparent resistivity Southern Line A,B={A},{B}")
plt.legend()
plt.xlabel("offset [m]")
plt.ylabel("rho [Ohmm]")
plt.savefig(os.path.join(PLOTDIR,"rho_south.png"))

plt.clf()
plt.scatter(X_north, np.array(data1_north)  , label = "i", s= 10 )
for kk in data_north:
    plt.scatter(X_north, np.array(data_north[kk]) , label = f"f{kk}", s= 10 )
plt.title(f"MN voltage Northern Line A,B={A},{B}")
plt.legend()
plt.xlabel("offset [m]")
plt.ylabel("voltage [V]")
plt.savefig(os.path.join(PLOTDIR,"data_north.png"))

plt.clf()
plt.scatter(X_south, np.array(data1_south) , label = "init", s= 10  )
for kk in data_south:
    plt.scatter(X_south, np.array(data_south[kk]) , label = f"f{kk}", s= 10 )
plt.title(f"MN voltage Southern Line A,B={A},{B}")
plt.legend()
plt.xlabel("offset [m]")
plt.ylabel("voltage [V]")
plt.savefig(os.path.join(PLOTDIR,"data_south.png"))

plt.clf()
for kk in data2_south:
    plt.scatter(X_south, np.array(data2_south[kk])*1000, label = f"f{kk}", s= 10 )
plt.title(f"MN voltage change Southern Line A,B={A},{B}")
plt.xlabel("offset [m]")
plt.ylabel("voltage [mV]")
plt.legend()
plt.savefig(os.path.join(PLOTDIR, "secondary_south.png"))

plt.clf()
for kk in data2_north:
    plt.scatter(X_north, np.array(data2_north[kk])*1000, label = f"f{kk}", s= 10 )
plt.title(f"MN voltage change Northern Line A,B={A},{B}")
plt.xlabel("offset [m]")
plt.ylabel("voltage [mV]")
plt.legend()
plt.savefig(os.path.join(PLOTDIR, "secondary_north.png"))

