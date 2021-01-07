# fingal - 3D Geophysical Inversion in python
A python package for geohpysical inversion in 3D for galvanic data based on the finite element method (FEM) supporting unstructured meshes. 
Data sets supported are electric potential data as well as electric field data for DC and IP surveys.  

It is based in the  [esys-escript](https://github.com/esys-escript/esys-escript.github.io) with python3 and supports parallel execution with domain decomposition for threading as well as MPI. 

# Installation 
For installation (sorry no setup.py yet) clone the `fingal` repository  

    mkdir fingal
    git clone https://github.com/LutzGross/fingal.git fingal

and add the path  to the `fingal` to your executable and python path:

    export PATH=${PATH}:${HOME}/fingal/bin
    export PYTHONPATH=${PYTHONPATH}:${HOME}/fingal/lib

Depending on your `esys-escript` installation you also have to add the path to 
the `esys-escript` installation to `PATH`, `LD_LIBRARY_PATH` and `PYTHONPATH`.

# Example

This is a simple example to demonstarte a  typical workflow from a file of electrodes (or stations) [electrodes.loc](examples/Example1/electrodes.loc) and data file of IP data [survey.csv](examples/Example1/survey.csv)  

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

The option `--optimize` rescales the reference conductivity `sigma0` set in the configuration file  prior to the inversion in an an attempt to reduce the inital misfit. `-d` switches on more output. The  
 
 
The switch `--xyx` activates the creation of a CSV file giving ccordinates and conductivity in core region set in the configuartion file. 

# List of Functions

# 


by @LutzGross
