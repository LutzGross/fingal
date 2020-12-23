# fingal - Geophysical inversion in python
A python package for geohpysical inversion for galvanic data based on the finite element method (FEM) supporting unstructured meshes. 
Data sets supported are electric potential data as well as electric field data for DC and IP surveys.  

It is based in the  [esys-escript](https://github.com/esys-escript/esys-escript.github.io) with python3 and supports parallel execution with domain decomposition for threading as well as MPI. 

# Installation 
For installation (sorry no setup.py yet) clone the `fingal` repository  

    mkdir fingal
    git clone https://github.com/LutzGross/fingal.git fingal

and add the path  to the `fingal` to your executable and python path:

    export PATH=${PATH}:${HOME}/fingal/bin
    export PYTHONPATH=${PYTHONPATH}:${HOME}/fingal/lib

Depending on your `esys-escript` installation you also have to ass the path to 
the `esys-escript` installation to `PATH`, `LD_LIBRARY_PATH` and `PYTHONPATH`.

# Example

First step is to create a configuration file `ex1.py`. The `mkConfig.py` script is helping to do this 
populating the file with some basic information:

    mkConfig.py --station electrodes.loc --data survey.csv ex1

The map of the station locations can be created using

    plotStations.py -i station.png ex1

    mkMesh 
# List of Functions

# 


by @LutzGross
