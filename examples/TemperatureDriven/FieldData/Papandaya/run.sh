export OMP_NUM_THREAD=5
cd t1; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..
cd t3; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..
cd t5; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..
cd t6; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..
cd t7; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..
cd t8; runTempInversion.py -d config 2>&1 | tee log.txt ; cd ..