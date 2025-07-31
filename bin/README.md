# Documentation 

## Scripts

(*not working* does mean that it is expected that the scripts could fixed easily.*)

### Meshing
- makeMesh.py  - create mesh from station data and topography (if needed)  
- mkMeshFromStations.py  - not sure this is still needed.
- mkMeshFromStationsWithTopo.py - not sure this is still needed.

- extractCoreSurfaceAsGeo.py - use this to create GMSH `.geo` file to add padding.

### Converter
- mph2fly.py - COMSOL  `.mph`-file to finley `.fly`-format
- vtu2fly.py - VTK `.vtu`-file to finley `.fly`-format
- gmsh2fly.py - GMSH `.msh`-file to finley `.fly`-format
- convertZZ.py - processess a ZZ-data file (not working)
- fly2silo.py - does exactly what the name says. 

### Inversions

- runSynthetic.py - creates a synethetic survey data set
- runSensitivityERT.py - sensitivity analysis for ERT surveys
- runERTinversion.py - ERT inversion
- runIPinversion.py - coupled IP inversion
- runIP2inversion.py - subsequent IP inversion from using ERT inversion results
- runTempInversion.py - temperature driven inversion (not working)
- runIPFieldinversion.py - FullWaver inversion (really not working) 

### Helpers
- mkWennerSurvey.py - makes a Wenner survey schedule file 
- plotStations.py - plot the stations
- mkSurveyFromINP.py -  INP file converter (not working) 
- mkSurveyFromGeosoftFile.py - GEOSOFT file converter (not working) 


## Configuration File

The configuration file is a *Python* script that sets the following variables:
 
### General information
 - `created` - some information on creator and creation date. 
 - `project` - project name

### Mesh
 
 - `meshfile` - name of the FEM mesh file of the domain for inversion, typically with extension `fly`.
 - `faces_tags` - list of tag names of the domain faces excluding the top surface. Typically, these are the faces of the domain that are 
orthogonal to the x-axis (left and right face) and the y-axis (front and bottom face) plus the bottom face. This is used by tools generating meshes.  
 - `surface_tags` - list of tag names of the (top) surfaces. This is used by tools generating meshes.   
 - `core_tags` - list of tag names defining the core region targeted by the inversion.
 - `padding_tags` - list of tag names defining the padding region added to the core region to eliminate boundary effects in the forward model and the inversion.

### Electrodes

 - `stationfile` -  name of the station (electrode) file in CSV format. Format is `i`, `x`, `y` or `i`, `x`, `y`, `z` where (`x`, `y`, `z`) is coordinate 
of the station with integer identifier `i`.
 - `stationdelimiter` - column seperator of the station file, e.g. `','`
 - `stationsFMT` - format string to convert station identifier `i` into the name tag for the coresponding FEM node, e.g. `'s%s'`.

### Data 

  - `dipoleInjections` -  Set `True` for ABMN and ABM style surveys.
  - `dipoleMeasurements` - Set `True` for AMN and ABMN style surveys. 
  - `datacolumns` - list specifying the values of the data file. The following tags are used:
    - `'R'` - resistance [V/A]
    - `'ETA'` - chargeability
    - `'ERR_R'` - error resistance [V/A]
    - `'ERR_ETA'` - error chargeability 
    - `'RELERR_R'` - rel. error resistance
    - `'RELERR_ETA'` - rel. error chargeability
    - `'E'` - electric field intensity per charge current [V/(Am)]
    - `'ERR_E'` - error electric field intensity per charge current [V/(Am)]
    - `'RELERR_E'` - rel. error electric field intensity per charge current
    
    If no error is set, a realtive error of 5% is used. For an IP survey with unknown errors, use `['R', 'ETA']`.
  - `datafile` - name of datafile in CSV format. Columns are the station identifiers of charging A and B (if applicable), the station identifiers for measurements M and N (if applicable) and the respective observations as declared by the `datacolumns` list.
  - `datadelimiter` - column seperator of the data file, e.g. `','`
  - `data_rtol` - drop tolerance for small values relative to the maximum value in the survey, e.g.  `1e-4`
  - `schedulefile` - name of survey schedule file in CSV format. The file is used when generating synthetic data. Columns are the station identifiers of charging stations A and B (if applicable) and the station identifiers for measurements M and N (if applicable).
  - `usesStationCoordinates` - if `True` station coordinates rather identifiers are used in the data file (not recommended to use), eg. `False`.

### Inversion


 - `sigma0_ref` - background conductivity in [S/m], e.g. `0.1`
 - `Mn_ref` - background normalised chargeability [S/m], e.g. `0.001`
   - `clip_property_function` - limitation value of the property function to avoid overflow in the expentential function 
   when calculating the conductivity, e.g. `10`.
 - `m_tolerance` - relative termination tolerance for the inversion on property function level, e.g. `1.e-4`
 - `g_tolerance` - relative termination tolerance for the inversion on cost function function level, e.g. `1.e-6`
 - `interpolation_order` - maximum interpolation order in , e.g. `3`. Must be less or equal three.
 - `imax` - maximum number of iterations, e.g. `400`
 - `truncation` - number of term in the BFGS schemes, e.g. `20`
 - `restart` - iteration count after which the BFGS is restarted (currently swiched off), e.g. `60`
 - `pde_tol` - tolerance for the PDE solver for the forward model, e.g. `1e-10`
 - `use_log_misfit_DC` - if `True` the log of the data is used in the DC misfit, e.g. `False`
 - `use_log_misfit_IP` - if `True` the log of the data is used in the IP misfit, , e.g. `False`
 - `regularization_weighting_DC_misfit` - weighting of the DC term in the misfit relative to the IP term, e.g. `0.5`. The factor is
is rescaled so the respective factors have the ratio `regularization_weighting_DC_misfit` and sum up to one.
 - `regularization_length_scale` - the length scale of the property function, e.g. `None`. If `None` infinity is used. 
This is only used for `H2` and `H2_0` regularization only.
 - `regularization_DC` - regularization type for the DC and IP inversion, e.g. `H1`. Allowed methods are `H1`, `H1_0`, `H2` and  `H2_0`:
   - `H1` : gradient of property function with fixed conductivity `sigma_ref`.  
   - `H1_0` : gradient of property function preserving conductivity mean `sigma_ref`
   - `H2` : Laplacian of property function with fixed conductivity `sigma_ref` on some boundaries 
   - `H2_0` : Laplacian of property function conductivity preserving mean `sigma_ref`: `H2_0`

   For the cases `H2` and `H2_0`,  setting the value `regularization_length_scale` has the effect that
  the gradient of the property function is added to the regularization with weighting factor `a**2` 
  for `a=1/config.regularization_length_scale`. 
   
 - `regularization_IP` - regularization type for IP2 inversion, e.g. `H1`. See `regularization_DC` for details.
 - `regularization_w1DC` - weighting factor for the regularization term for the log of conductivity in the cost function, e.g. `1e-1`
 - `regularization_w1IP` -  weighting factor for the regularization term for the log of normalized chargeabilty in the cost function, e.g. `1e-1`
 - `regularization_theta` - weighting factor for the cross-gradient term in the cost function, e.g. `0.`
 - `use_robin_condition_in_model` - if `True`, Robin boundary condition in the forward model is applied, e.g. `False`. This allows for using smaller domains but requires resolving PDEs for each injection. In general, this model 
is computationally more expensive in comparison of using a larger domain with zero potential on the boundray. Not all inversions support this model (yet).
 - `use_L1Norm` - if `True` L1-norm is used in the `H1`-regularization, e.g. `False`. **It needs more testing.**
 - `epsilon_L1Norm` - cut-off tolerance for the L1-norm in the regularization, e.g. `1e-8`
 - `fixed_region_tags` - list of tags of domain  faces where physical properties are fixed, using the tags introduced by
`faces_tags` and `surfaces_tags`, e.g. `[]`. If `None` or `[]` and boundary conditions are required, then these are the faces of the domain that are 
orthogonal to the x-axis (left and right face) and the y-axis (front and bottom face) plus the bottom face. Notice that in this case it is assumed that these faces are plain.
 - `fix_top` - If `True` and `fixed_region_tags` is `None`, in addition the top surface is added,e.g. `False`. 
### Output handeling:

 - `outfile`- name of the result file, e.g.`'sigma'`. Depending on the requested output format an appropriate extension is added.
 - `sigma0_dump` - name of dump file for DC conductivity, e.g. `'sigma0dump'`. If set, ERT inversion is creating this file to provide the low frequency conductivity for a seperate IP2 inversion.

### True Properties

The synthetic data generator expects a function of the name `true_properties` that delivers the true electric conductivity and the true normalised 
chargeability at a given domain. This is an example:

    def true_properties(domain):
        from esys.escript import Scalar, Function
        sigma0_true=Scalar(sigma0_ref , Function(domain))
        sigma0_true.setTaggedValue('anomaly', sigma0_ref * 10)
        Mn_true=Scalar(Mn_ref, Function(domain))
        Mn_true.setTaggedValue('anomaly',  sigma0_ref * 10 * 0.5 )
        return sigma0_true, Mn_true


