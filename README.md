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

# Examples

- As simple example for a 2D ERT inversion: [Example1](examples/Example1/README.md)
- Synthetic FullWaver data inversion: [SyntheticFW](examples/SyntheticFW/README.md)
- Volcano with FullWaver data: [VolcanoFW](examples/VolcanoFW/README.md)

# List of Functions

# 


by @LutzGross
