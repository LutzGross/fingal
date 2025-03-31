# Mine

## Mesh generation

Generate mesh form `mine3D.geo:`

    ./mkMesh

This read `stations.csv`, renerates a mesh file `mine.msh` using `gmsh` and 
converts the mesh file `mine.fly` in `fly` format. 

## Synthetics

Create the sensitivity and resolution loss map:

    runSensitivityERT.py config
