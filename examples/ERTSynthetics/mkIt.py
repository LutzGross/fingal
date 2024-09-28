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
NumberOfLines = 1
distanceLines = 0
with open(GEO_FILE, 'r' ) as f:
    for line in f:
        if line.startswith('numElectrodes'):
           NumberOfElectrodes =int( line.split("=")[1][:-2] )
        if line.startswith('distanceElectrodes'):
            distanceElectrodes= float(line.split("=")[1][:-2] )
        if line.startswith('numLines'):
           NumberOfLines =int( line.split("=")[1][:-2] )
        if line.startswith('distanceLines'):
            distanceLines= float(line.split("=")[1][:-2] )

assert not distanceElectrodes is None
assert not NumberOfElectrodes is None
LengthOfLine = distanceElectrodes * (NumberOfElectrodes-1)
WidthOfLines = distanceLines * (NumberOfLines-1)

print(f"searching file {GEO_FILE}.")
print(f"number of electrodes per line = {NumberOfElectrodes}")
print(f"electrode spacing in line = {distanceElectrodes}")
print(f"length of line = {LengthOfLine}")
print(f"number of lines = {NumberOfLines}")
print(f"spacing of line = {distanceLines}")
print(f"width of lines = {WidthOfLines}")
# ------


# ... create stations
k1=0
stations={}
for i in range(NumberOfLines):
    for j in range(NumberOfElectrodes):
        xpos=-LengthOfLine/2+distanceElectrodes*j
        ypos=-WidthOfLines/2+distanceLines*i
        ipos = 1000 * (i+1) + (j+1)
        stations[ipos+k1]=(xpos, ypos, 0.0)
f=open(config.stationfile,'w')
for s in stations:
    f.write("%d, %s, %s, %s\n"%(s, stations[s][0], stations[s][1], stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(len(stations), config.stationfile ))

#... create schedule:
schedule=[]
for i in range(NumberOfLines):
    schedule+=makeWennerArray(numElectrodes=NumberOfElectrodes, id0=k1+1000 * (i+1) + 1)

f=open(config.schedulefile,'w')
for A, B, M, N in schedule:
    f.write("%s, %s, %s, %s\n"%(A, B,M, N))
f.close()
print("%s scheduled observations written to file %s"%(len(schedule),config.schedulefile))
