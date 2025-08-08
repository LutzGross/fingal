# Synthetic IP Inversion with DC and IP solved separately 

### Generate Data Set
The first step is to create a synthetic data set. The [gmsh](https://gmsh.info/)  
geometry file defines core domain with some padding and subdomains staged as `anomaly_right`
and `anomaly_left` that define electric conductivity anomalies
for the synthetic data generation. 
Edit ths `geo`-file [`with_anomaly.geo`](./with_anomaly.geo) to change the shape and 
location of the anomaly and the survey that is set as equidistant, parallel lines of equidistant
electrodes.

Run 

    gmsh -3 -optimize_netgen -o mesh_synthetic.msh with_anomaly.geo

to generate the mesh for the synthetic data generation. 
To make the station and schedule file (names are set in [`config.py`](./config.py)) run

    python3 mkIt.py

The script inspects the `with_anomaly.geo` file grab the number of electrodes, number of lines
and their spacings. The schedule is following a Wenner set-up for each line. 


### Run Inversion

Plot the stations to file `plot_stations.png`:

    plotStations.py --image plot_stations --debug config


Create the data file with 1% noise:

    runSynthetic.py --noise 1 --silo syntheticmesh mesh_synthetic.msh config


This is a simple way to generate a mesh from the electrode's positions:

    mkMeshFromStations.py --coremeshfactor 0.5 config

The mesh is written to the `config.meshfile` in the *esys.finley* `fly` format.

First step in the inversion is to run an ERT inversion:

    runERTinversion.py --vtk -d config

We refer to [SyntheticsERT](../SyntheticsERT/README.md) for details.
The inversion - if successful - will create a dump file `dump_sigma0` as been set
in `config.sigma0_dump` which provides the input for the low frequency
conductivity `sigma` for the IP inversion targeting the normalised chargeability 
`Mn`. 


    runIP2inversion.py --vtk -d config

 


by @LutzGross
