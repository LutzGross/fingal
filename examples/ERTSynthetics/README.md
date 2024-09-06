# Synthetic ERT Inversion Example

### Generate Data Set
The first step is to create a synthetic data set. The [gmsh](https://gmsh.info/)  
geometry file defines core domain with some padding and subdomain stagged as `anomaly` 
that is used to
define an electric conductivity anomaly for the synthetic data generation. 
Edit this file to change the shape and location of the anomaly. 

Run 

    gmsh -3 -optimize_netgen -o mesh_synthetic.msh with_anomaly.geo

to generate the mesh for the synthetic data generation. 
To make the station and schedule file (names are set in [`config.py`](./config.py)) run

    python3 mkIt.py

The script inspects the `with_anomaly.geo` file grab the number of electrodes and the spacing. 
The schedule is following a Wenner set-up. This can be changed if needed.

Create the data file with 5% noise:

    runSynthetic.py --noise 5 --silo syntheticmesh mesh_synthetic.msh config



=======================

    cd examples/Example1
    mkConfig.py --station electrodes.loc --data survey.csv --obs R,ERR_R,ETA ex1

To create map of the station locations can be created using

    plotStations.py -i station.png ex1

The [gmsh](https://gmsh.info/) mesh generator is used to create an esys-escript fly file of a 3D finite element mesh: 

    mkMesh.py --geo mymesh ex1 

The file name is `ex1.fly` as defined in the configuration file `ex1.py`. You can use gmsh to inspect the  geometry file `mymesh.geo` and mesh `mymesh.msh`. 
The geometry is defined as the area spanned by the electrodes plus a little extra area. This core region which is vertically extended is meshed with a finer mesh where extra refinement near the position of the electrodes is added. Around the core region additonal padding is introduced to allow for a smooth deacy of the electrial potentials towards the boundary. This extra area is meshed with a courser mesh. `mkMesh.py` allows to configure this geometrical set-up, see `mkMesh.py -h` for details.

To run the inversion based on the configuration file `ex1.py` use: 

    runERTinversion.py --optimize   --xyz --vtk -d ex1

The option `--optimize` rescales the reference conductivity `sigma0` set in the configuration file  prior to the inversion in an an attempt to reduce the inital misfit. `-d` switches on more output. With the switch `--vtk` the file `sigma.vtk` in the [VTK](https://vtk.org/) file format is created where `sigma` is taken from the `output` variable set in `ex1.py`  (by default the `silo` format is used which is more compact but less portable). You can use 3D visualization packages such as

- [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit)
- [paraview](https://www.paraview.org/)
- [mayavi](https://docs.enthought.com/mayavi/mayavi/)
 
The result should look this: 
<p align="center">
  <img src="result.png" width="600" title="VisIt Visualization">
</p>

The switch `--xyx` activates the creation of a CSV file giving ccordinates and conductivity in core region set in the configuartion file via the `core` variable. The name of the created file is `sigma.csv` where again `sigma` is taken from the `output` variable set in `ex1.py`. You can read and plot this file for instance using a 3D scatter plot in `matplotlib`:

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d
    data= np.loadtxt('sigma.csv',skiprows=1, delimiter=',')
    volume=data[data[:,3]>0.04]
    xyz=volume[:,:3]
    sigma=volume[:,3]
    ax = plt.axes(projection='3d')
    ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], c = sigma, s=0.01)
    plt.show()


by @LutzGross
