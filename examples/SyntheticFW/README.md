# Synthetic Field Intensity Inversion (aka [FullWaver](http://www.iris-instruments.com/v-fullwaver.html))

This is synthetic case with three cubic anomalies where
one cube defines anomalies in conductivity and charageability,
one cube defines an anomaliy in conductivity and
one cube defines an anomaliy in charageability. The doamin is defines in the Gmsh file `domain.geo` which can be visualized in gmsh: 

    cd examples/SuntheticFW
    gmsh domain.geo    
    
and show this geometry setup:
<p>
    <img src="domain.png" width="600" title="geometry of the domain">
</p>
The core part shows the three anomalies and the position of charging electrodes and recording stations. In the core region a finer mesh is used while
in the padding region which is added to eliminate boundary effects a coarser mesh is used. We will create a synthetic data set containing 
electric field intensity and (modified) chargeability [Missing reference](XX). 

The first step is to extract the location of the charging electrodes and recording stations from the `domain.geo` to create the file of station locations `stations.csv` as defined in the configuration file [`config.py`](config.py) and to generate a schedule file `synth_data.csv`  (again the name is defined in the configuration file) which will later be used to run the virtual survey: 

    ./mkSchedule.py -d 10 domain.geo config

In the survey the electrodes at an offset in the East, North, West and South are paired with in this case 10 (specified by the option `-d 10`) randomly choosen other electrodes setting in this case 4 x 10 charging experiments. The station spacing is recorded as about 300 (if you have not changed the geometry file). 

The lacation of the stations and electrodes can be plotted to the file `station_positions.png` uisng  

    plotStations.py -i station_positions.png config

where green and blue dots refer to recording stations and charging electrode respectively:     
<p>
    <img src="station_positions.png" width="600" title="Position of measurement stations and chargong electrodes">
</p>
A 3D mesh is generated from `domain.geo` using gmsh with the mesh written to `synth.msh` in the GMSH file format (Make sure that the format 2.2 is used):

    gmsh -3 -o synth.msh domain.geo
    
Notice the the cubic anomlies are tagged with "Anomaly1", "Anomaly2" and "Anomaly3". The core region excluding the anomalies is tagged "InnerBox"
while the padding is the tagged "OuterBox". The mesh will have about 135000 nodes. The resolution can be changed by editing `domain.geo`.

The gmsh file ' synth.msh' is converted into an [esys-escript](https://github.com/esys-escript/esys-escript.github.io) mesh file mainly for reasons of performance when reading the mesh under MPI: 

    gmsh2fly.py --silo mesh synth.msh config

The output file name in the `fly` format is specifeid in the configuration file [`config.py`](config.py). The converter aslo generates `mesh.silo` to visualize mesh as it is  
<p>
    <img src="mesh.png" width="600" title="Position of measurement stations and chargong electrodes">
</p>

....

This is a simple example to demonstarte a  typical workflow from a file of electrodes (or stations) [electrodes.loc](examples/Example1/electrodes.loc) and data file of IP data [survey.csv](./survey.csv)  

First step is to create a configuration file `ex1.py`. The `mkConfig.py` script is helping to do this 
populating the file with some basic information:

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


by @LutzGross
