# .. Instructions ...
#
# ...  create station file and schedule:
#
#    ./mkIt.py
#
# ... generate mesh fly mesh
#
#   mkMeshFromStations.py config
#
# ... generate synthetic data
#
#   runSynthetic.py -d --silo survey config
import importlib, sys, os
from fingal  import makeWennerArray

CONFIG_FILE='config'
# ----
LengthOfLine = 80
NumberOfElectrodes = 64
ElectrodeDistance = LengthOfLine/(NumberOfElectrodes-1)
# ------

config = importlib.import_module(CONFIG_FILE)
print("configuration "+CONFIG_FILE+" imported.")

# ... create stations
k1=1
stations={}
for j in range(NumberOfElectrodes):
        xpos=ElectrodeDistance*j
        ypos=0.
        stations[10*j+k1]=(xpos, ypos, 0.0)

f=open(config.stationfile,'w')
for s in stations:
    f.write("%d, %s, %s, %s\n"%(s, stations[s][0], stations[s][1], stations[s][2] ))
f.close()
print("%s stations  written to file %s"%(NumberOfElectrodes, config.stationfile ))

#... create schedule:
schedule=makeWennerArray(numElectrodes=NumberOfElectrodes, id0=0)

f=open(config.schedulefile,'w')
for A, B, M, N in schedule:
    f.write("%s, %s, %s, %s\n"%(10*A+k1, 10*B+k1,10*M+k1, 10*N+k1))
f.close()
print("%s schedules observations written to file %s"%(len(schedule),config.schedulefile))
