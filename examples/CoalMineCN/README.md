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

It requires `../mine.fly` and `../stations.csv` and creates `data.csv` from `schedule.csv`.
To run the inversion.

    runERTinversion.py --nooptimize -d config


## Field Data

    cd FieldData


It uses `../mine.fly`, `../stations.csv` and  `data_field.csv`.


To run the inversion.

    runERTinversion.py --nooptimize -d config

