# Example: ERT Field Data 
 
This is an example of the inversion of an ERT data set collected on Heron Island of the 
Queendland Central Coast. It is a single Schlumberger line of 64 electrodes 
on a assumed flat surface. 

To run the inversion two data set are required:

- **electrode/station locations** in `csv` format with each rom giving station identifier 
and 3D station coordinates. (note: for this case it is assumed that all z/x_2 coordinates are identical). 
Ideally the line should be aligned with x or y axis.
- **data file** also in `csv` format. 
