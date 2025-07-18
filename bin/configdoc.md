# Configuration File

The configuration file is a *Python* script that sets the following variables:
 
### General information
 - `created` - some information on creator and creation date. 
 - `project` - project name

### Mesh
 
 - `meshfile` - name of the FEM mesh file of the domain for inversion, typically with extension `fly`.
 - `faces_tags` - list of tag names of the faces at which properties are fixed. Typically, these are the faces of the domain that are 
orthogonal to the x-axis (left and right face) and the y-axis (front and bottom face) plus the bottom face.     
 - `surface_tags` - list of tag names of the faces at which properties are updated. Typically, this is the top surface.   
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


 - `sigma0_ref` - , e.g. `0.1`
 - `Mn_ref` - , e.g. `0.001`
 - `fixed_region_tags` - , e.g. `[]`
 - `fix_top` - , e.g. `False`
 - `clip_property_function` - , e.g. `10`
 - `m_tolerance` - , e.g. `1.e-4`
 - `g_tolerance` - , e.g. `1.e-6`
 - `interpolation_order` - , e.g. `3`
 - `imax` - , e.g. `400`
 - `truncation` - , e.g. `20`
 - `restart` - , e.g. `60`
 - `pde_tol` - , e.g. `1e-10`

 - `use_log_misfit_DC` - , e.g. `False`
 - `use_log_misfit_IP` - , e.g. `False`

 - `regularization_weighting_DC_misfit` - , e.g. `1`
 - `regularization_DC` - , e.g. `'H1'` # in ['H1', "H1_0", 'H2',  "H2_0"]
 - `regularization_IP`
 - `regularization_length_scale` - , e.g. `None` # only used for "H2" and "H2_0" regularization
 - `regularization_w1DC` - , e.g. `1e-1`
 - `regularization_w1IP` - , e.g. `1e-1`
 - `regularization_theta` - , e.g. `0.`
 - `use_robin_condition_in_model` - , e.g. `False`
 - `use_L1Norm` - , e.g. `False`
 - `epsilon_L1Norm` - , e.g. `1e-8`
 
### Output handeling:

 - `outfile`=`'sigma'`

 - #restartfile = 'restart'

### True Properties


    def true_properties(domain):
        from esys.escript import Scalar, Function
        sigma0_true=Scalar(sigma0_ref , Function(domain))
        sigma0_true.setTaggedValue('anomaly', sigma0_ref * 10)
        Mn_true=Scalar(Mn_ref, Function(domain))
        Mn_true.setTaggedValue('anomaly',  sigma0_ref * 10 * 0.5 )
        return sigma0_true, Mn_true


