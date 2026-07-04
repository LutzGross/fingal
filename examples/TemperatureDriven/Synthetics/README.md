# Synthetic Temperature-driven Inversion

Edit [sonfig.py](config.py) 

Generate Wanner schedule and stations:

    mkWennerSurvey.py --numElectrodesPerLine 64 --numLines 1 --lineLength 320 --lineDistance 10 config

Generate Data:

     ./mkData.py --noise 0 --silofile truemodel config 

Make inversion mesh:

    makeMesh.py --fly mesh.fly stations.csv 

Run inversion:

      




