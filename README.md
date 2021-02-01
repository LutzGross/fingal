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
    export PYTHONPATH=${PYTHONPATH}:${HOME}/fingal/bin

Depending on your `esys-escript` installation you also have to add the path to 
the `esys-escript` installation to `PATH`, `LD_LIBRARY_PATH` and `PYTHONPATH`.

# Examples

- As simple example for a 2D ERT inversion: [Example1](examples/Example1/README.md)
- Synthetic FullWaver data inversion: [SyntheticFW](examples/SyntheticFW/README.md)
- Volcano with FullWaver data: [VolcanoFW](examples/VolcanoFW/README.md)

# List of Command Line Functions

This is a list of some of the command line function provided. Use the command line option `-h` to see the options.

- `mkConfig.py` : creates an initial configuartion file that need to be adjudted to the specific case.
 - `mkMesh.py`: uses a station file to set up a gmsh geo file to generate finite element mesh
- `runSynthetic.py`: run a synthetic survey using a given schedule. 
- `convertZZ.py`: converts a ZZ survey file into a station and data file.
- `fly2silo.py`: creates a silo 3D visualization file from fly mesh file.
- `gmsh2fly.py`: converts a gmsh msh file to fly mesh file.
- `mph2fly.py`: converts a COMSOL mph file to fly mesh file.
- `plotStations.py`: plot the survey stations in plane view. 
- `runERTinversion.py`: runs ERT (an IP) inversion using potential fields.
- `runIPFieldinversion.py`: runs ERT (an IP) inversion using electric field data (aka FullWaver).

by @LutzGross
