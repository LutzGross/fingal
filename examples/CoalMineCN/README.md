# Mine

In this application two ERT surveys are run on gallery walls in a mine, see [`mine3D.geo`](./mine3D.geo) for the geometry.
There are two 



## Mesh generation

Generate mesh form `mine3D.geo:`

    ./mkMesh

This read `stations.csv`, renerates a mesh file `mine.msh` using `gmsh` and 
converts the mesh file `mine.fly` in `fly` format. 

## Sensitivity

Create the sensitivity and resolution map:

    runSensitivityERT.py config

## Synthetics

Create synthetic survey run

    cd Synthetics
    ./mkData.py

creat


