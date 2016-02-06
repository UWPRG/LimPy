Instructions for how to use Langevin Integrator:
LangevinFunctions.py contains the functions for the actual integrator
along with the functions defining the potential energy/force calculations and
functions defining what is a rare event for infrequent sampling.

LangevinScript is what should be run to implement everything. It sets up a
GUI for inputs from the user. Also, inputs can be read in from a CSV file if
the box is checked and file name provided.

As of now the potential energy surfaces are hard coded in (as these were the
only two I had). The functions can implement MD, Metadynamics, Well-Tempered
Metadynamics, and infrequent metadynamics.

The script will create two output files (and a series of pngs of plots if
flagged). One output file will have the data for the infrequent Metadynamics
ipython notebook. This output file name is specified by the user. The other
output file is the same as the other output file but has 'info' added to the
name. This file has the input parameters used, as well as RMSD values if the
test was run for a converged PES. 
