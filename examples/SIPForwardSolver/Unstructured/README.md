These are test cases for unstructured meshes.

First  generate a series of meshes in `fly`-format :

    mkdir meshes
    mkManyMeshes.py

It is a source `"S1"` over a three anomalies 
`"Anomaly1"`, `"Anomaly2"` and `"Anomaly3"`.

To run the test example with conductivity contrast of 10 and a ratio of
0.05 of imaginary over real part of conductivity in the anomalies:

    ./unstructest.py  --contrast 10 --ratio 0.05  meshes/dom16.fly 

Add `--test` to compare the iterative solver solution with the direct solver solution. This should be used for small meshes only.

Create a script `run.sh` to run tests with different meshes, values of `ratio` and values of `contrast`:

    ./mkRuns.py
    mkdir logs

Log runs are written into directory logs. 
