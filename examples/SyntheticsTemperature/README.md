# Synthetic Temperature-driven Inversion

Edit [sonfig.py](config.py) 

Generate Wanner schedule and stations:

    mkWennerSurvey.py --numElectrodesPerLine 32 --numLines 5 --lineLength 320 --lineDistance 10 config

Generate Data:

    mkData.py --noise 0 --coredepth 45 --extracore 20 --padding 150 --coremeshfactor 1. --stationmeshfactor 0.3 --paddingmesh 20 --silofile truemodel config    

Make inversion mesh:

    mkMeshFromStations.py --coredepth 45 --extracore 20 --padding 150 --coremeshfactor 1. --stationmeshfactor 0.3 --paddingmesh 20 --silofile mesh config



