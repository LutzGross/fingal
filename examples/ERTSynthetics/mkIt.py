# .. Instructions ...
#
# ...  create station file and schedule:
#
#    ./mkIt.py
#
# ...
#
import importlib, sys, os
from fingal  import makeWennerArray


CONFIG_FILE = 'config'
GEO_FILE = 'with_anomaly.geo'

config = importlib.import_module(CONFIG_FILE)
print("configuration "+CONFIG_FILE+" imported.")
# find number of electrodes and Distance of electrodes in GEO file:
# ----
NumberOfElectrodes = None
distanceElectrodes = None
with open(GEO_FILE, 'r' ) as f:
    for line in f:
        if line.startswith('numElectrodes'):
           NumberOfElectrodes =int( line.split("=")[1][:-2] )
        if line.startswith('distanceElectrodes'):
            distanceElectrodes= float(line.split("=")[1][:-2] )

assert not distanceElectrodes is None
assert not NumberOfElectrodes is None
LengthOfLine = distanceElectrodes * (NumberOfElectrodes-1)

print(f"searching file {GEO_FILE}.")
print(f"number of electrodes = {NumberOfElectrodes}")
print(f"electrode spacing = {distanceElectrodes}")
print(f"length of line = {LengthOfLine}")
# ------


# ... create stations
k1=1
stations={}
for j in range(NumberOfElectrodes):
        xpos=-LengthOfLine/2+distanceElectrodes*j
        ypos=0.
        stations[j+k1]=(xpos, ypos, 0.0)

f=open(config.stationfile,'w')
for s in stations:
    f.write("%d, %s, %s, %s\n"%(s, stations[s][0], stations[s][1], stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(NumberOfElectrodes, config.stationfile ))

#... create schedule:
schedule=makeWennerArray(numElectrodes=NumberOfElectrodes, id0=k1)

f=open(config.schedulefile,'w')
for A, B, M, N in schedule:
    f.write("%s, %s, %s, %s\n"%(A, B,M, N))
f.close()
print("%s scheduled observations written to file %s"%(len(schedule),config.schedulefile))
