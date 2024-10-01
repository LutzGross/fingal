# Example: ERT Field Data 
 
This is an example of the inversion of an ERT data set collected on Heron Island of the 
Queendland Central Coast. It is a single Schlumberger line of 64 electrodes 
on a assumed flat surface. 

To run the inversion two data sets are required:

- **electrode/station locations** in `csv` format with each rom giving station identifier 
and 3D station coordinates. (note: for this case it is assumed that all z/x_2 coordinates are identical). 
Ideally the line should be aligned with x or y axis.
- **data file** also in `csv` format. Each row gives A, B, M, N (using identifyers) plus
the respective measurements. Order of data columens are specified be list `datacolumns` (here just `R`) 
with the following meanings (if an error is not given 5% is used):

  - 'R' - resistance [V/A]
  - 'E' - electric field intensity per charge current [V/(Am)]
  - 'ETA' - chargeability
  - 'ERR_R' - error resistance [V/A]
  - 'ERR_E' - error electric field intensity per charge current [V/(Am)]
  - 'ERR_ETA' - error chargeability 
  - 'RELERR_R' - rel. error resistance
  - 'RELERR_E' - rel. error electric field intensity per charge current
  - 'RELERR_ETA' - rel. error chargeability  [1]

In addition, a configuration file is requested. It is a *Python* script where an 
initial version can be generated from the data and station file using 
the `mkConfig.py` script.   

## Workflow

First step is to generate a mesh from the station file specified in the configuration file
[`config.py`](/config.py):

    mkMeshFromStations.py config

This generates the mesh file `mesh.fly`  in *esys.finley* mesh file (file name is 
set by `meshfile` variable in the configuration file). The mesh written isd also
written to `mesh.silo` for visualization.

Plot of the stations? Use:

    plotStations.py --image plot_stations --debug config

Run the inversion

    runERTinversion.py --vtk -d config  

See [ERTSynthetics/README.md](../ERTSynthetics/README.md) for more details.

by @LutzGross