# Synthetic ERT Inversion Example

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

to generate the mesh for the synthetic data generation. To make the station and schedule file (names are set in [`config.py`](./config.py)) run

    python3 ./mkIt.py

The script inspects the `with_anomaly.geo` file grab the number of electrodes, nuber of lines
and their spacings. The schedule is following a Wenner set-up for each line. 


### Run Inversion

Plot the stations to file `plot_stations.png`:

    plotStations.py --image plot_stations --debug config


Create the data file with 1% noise:

    runSynthetic.py --noise 1 --silo syntheticmesh mesh_synthetic.msh config


This is a simple way to generate a mesh from the electrode's positions:

    mkMeshFromStations.py --coremeshfactor 0.5 config

The mesh is written to the `config.meshfile` in the *esys.finley* `fly` format.

To run the inversion based on the configuration file `config.py` use: 

    runERTinversion.py --vtk -d config

The reference conductivity `sigma_ref` set in the configuration file 
is reset prior to the inversion in an attempt to reduce the initial misfit.  
(use line command `--nooptimize` to switch off this function).
`-d` switches on more output. With the switch `--vtk` the file `sigma.vtk` in 
the [VTK](https://vtk.org/) file format is created where `sigma` is taken from the `output` 
variable set in configuration file
(by default the `silo` format is used which is more compact 
but less portable). 

Misfit can be squares of data residual (weighted by error) (aka as chi^2) or square of 
data logarithm (`config.use_log_misfit_ERT` = `True`). 

Various regularization approaches can be used: 

   - Gradient of property function with fixed conductivity `sigma_ref` on some boundaries: `H1`
   - Gradient of property function preserving conductivity mean `sigma_ref`: `H1_0`
   - Laplacian of property function with fixed conductivity `sigma_ref` on some boundaries: `H2`, 
   - Laplacian of property function conductivity preserving mean `sigma_ref`: `H2_0`

For the cases `H2` and `H2_0` one can also set `config.regularization_length_scale` with the effect that
the gradient of the property function is added to the regularization with weighting factor `a**2` 
for `a=1/config.regularization_length_scale`.  


## Sensitivity and Resolution

Files of a sensitivity density and resolution loss map can be created by

    runSensitivityERT.py --file sensitivity --obs 1001,1032,1008,1024 config

where `--file`  is the file name for the _silo_ output file (use option `--vtk`
for VTK file format). The sensitivity density for a specific experiment
The mesh is read from the mesh file given in the configuration file.
Alternatively use the `--mesh` option (in `msh` or `fly` format). It is
assumed that conductivity is constant. When `--truesigma` is set using 
the function `true_properties` in the configuration file if set.

For the sensitivity density _S_ at a point and a change _dsigma_ 
to conductivity volume covering a volume _V_ 
the mean apparent conductivity over the survey is changing by 
_S * V * dsigma_ . The value of the resolution loss at a point 
specifies the edge length of a conductivity around that point that 
would create the same change to the mean apparent conductivity 
as at the location of highest sensitivity density, typically ot/near 
the electrodes. 


by @LutzGross
