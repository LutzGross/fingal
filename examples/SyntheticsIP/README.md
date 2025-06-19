# Synthetic IP Inversion Example

### Generate Data Set
The first step is to create a synthetic data set. The [gmsh](https://gmsh.info/)  
geometry file defines core domain with some padding and subdomains stagged as `anomaly_right`
and `anomaly_left` that define electric conductivity anomalies
for the synthetic data generation. 
Edit ths `geo`-file ['with_anomaly.geo'](./with_anomaly.geo) to change the shape and 
location of the anomaly and the survey that is set as equidistant, parallel lines of equidistant
electrodes.

Run 

    gmsh -3 -optimize_netgen -o mesh_synthetic.msh with_anomaly.geo

to generate the mesh for the synthetic data generation. 
To make the station and schedule file (names are set in [`config.py`](./config.py)) run

    python3 mkIt.py

The script inspects the `with_anomaly.geo` file grab the number of electrodes, nuber of lines
and their spacings. The schedule is following a Wenner set-up for each line. 


### Run Inversion

Plot the stations to file `plot_stations.png`:

    plotStations.py --image plot_stations --debug config


Create the data file with 1% noise:

    runSynthetic.py --noise 1 --silo syntheticmesh mesh_synthetic.msh config


This a simple way to generate a mesh from the electrode's positions:

    mkMeshFromStations.py --coremeshfactor 0.5 config

The mesh is written to the 'config.meshfile' in the *esys.finley* `fly` format.

To run the inversion based on the configuration file `config.py` use: 

    runIPinversion.py --vtk -d config

The reference conductivity `sigma_ref` and normalised chargeability are set in
the configuration file 
is rescaled prior to the inversion in an attempt to reduce the initial misfit
(use '--nooptimize' to switch off rescaling).
`-d` switches on more output. With the switch `--vtk` the file `sigma.vtk` in 
the [VTK](https://vtk.org/) file format is created where `sigma` is taken from the `output` 
variable set in configuration file
(by default the `silo` format is used which is more compact 
but less portable). 

Misfit can be squares of data residual (weighted by error) (aka as chi^2) or square of 
data logarithm (set `use_log_misfit_ERT = True` and `use_log_misfit_IP = True`). 

Various regularization approaches can be used: 

   - Gradient of property function: `'H1'`
   - Gradient of property function with zero-mean constraint: `'H1_0'`
   - Laplacian of property function: `'H2'`
   - Gradient of property function with zero-mean constraint: `'H2_0'`
 
REVISE
For H2-regularization uses the addition regularization parameter
`regularization_length_scale`. If set (not equal to `None`) then the value defines the 
length scale (in a sense of the exponential variogram). If not set (=`None`)
length scale infinity is used. 


by @LutzGross
