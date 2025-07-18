# Example: ERT Field Data 
 
This is an example of the inversion of an ERT data set collected on Heron Island of the 
Queensland Central Coast. It is a single Schlumberger line of 64 electrodes 
on an assumed flat surface. 

To run the inversion, two data sets are required:

- **electrode/station locations** in `csv` format with each row giving the station identifier 
and 3D station coordinates. (Note: for this case, it is assumed that all z/x_2 coordinates are identical.) 
Ideally, the line should be aligned with x or the y-axis.
- **data file** also in `csv` format. Each row gives A, B, M, N (using identifiers) plus the respective measurements. The order of data columns is specified by the list `datacolumns` (here just `R`), see [configuration files](../../bin/configdoc.md). 

In addition, a configuration file is requested. It is a *Python* script where an 
initial version can be generated from the data and station file using 
the `mkConfig.py` script.   

## Workflow

The first step is to generate a mesh from the station file specified in the configuration file
[`config.py`](/config.py):

    mkMeshFromStations.py config

This generates the mesh file `mesh.fly`  in *esys.finley* mesh file (file name is 
set by `meshfile` variable in the configuration file). The mesh written isd also
written to `mesh.silo` for visualisation. 

**Note**: The mesh is generated using [gmsh](https://gmsh.info/). The geometry file and *gmsh* mesh file can be inspected, `tmp.geo` and `tmp.msh` respectively.

Plot of the stations? Use:

    plotStations.py --image plot_stations --debug config

Run the inversion

    runERTinversion.py --vtk -d config  

See [ERTSynthetics/README.md](../SyntheticsERT/README.md) for more details.

by @LutzGross
