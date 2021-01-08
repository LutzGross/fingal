# Synthetic Field Intensity Inversion (aka FullWaver)

The 

    cd examples/SuntheticFW
    gmsh domain.geo    
<p>
    <img src="domain.png" width="600" title="geometry of the domain">
</p>

Extract stations form the geo file and create a survey schedule file:

    ./mkSchedule.py -d 10 domain.geo config

Plot the stations:

    plotStations.py -i station_positions.png config
    
<p>
    <img src="station_positions.png" width="600" title="Position of measurement stations and chargong electrodes">
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
