# Langevin_Package
----
Langevin_Package is a Langevin Integrator capable of carrying out Molecular Dynamics
(MD), Metadynamics, Well-Tempered Metadynamics, and Infrequent Metadynamics.
The integrator is capable of operating on, and visualizing, 1-D and 2-D potentials.
This package can be operated to run from the command line (`run_langevin.py`) or through a GUI (`run_langevin_GUI.py`). A version of this code has also been created with MPI4Py so it may be implemented in parallel on multiple nodes (`run_langevin_parallel.py`)

----
## Quick Start:
----
### Inputs:
The inputs specifying the system, method, parameters, and visualization options can be received via an inputs file or through a user interface.

The inputs files contains all of the specified inputs for the system of interest in the format below. Note information following #'s are comments indicating other options for those inputs.

```
Parameter,Input
Dimension,1-D Potential#2-D Potential
Method,Well-Tempered Metadynamics#Well-Tempered Metadynamics, MD, Infrequent WT MetaD, Metadynamics
Potential_Function,cosine_potential
Plotting,True
Steps,100000
Step size,0.01
Temperature,300
Mass,1
X0,3.14159
Xmin,-0.31
Xmax,6.5
Xincrement,0.01
Y0,0.5
Ymin,-2
Ymax,2
Yincrement,0.01
Gaussian Height,0.25
Gaussian Width,0.07
Deposition Frequency,250
Well Temperature,10000
Trials,100
Data Filename,output
Gamma,5
Plot Freq,1000
Make Movie,True
Kb,0.001987#kcal/mol
```
**All lines must be included for a properly formatted input file, however a line is not applicable to a system (ie Y parameters for a 1-D system), then just ignore those contents are change to values to 0.**

Each of the definitions can be defined in the key below.

----
 ![alt text](https://github.com/UWPRG/Chris_Scripts/blob/master/Langevin_Package/Input%20Definitions.png)
----

----
#### Parallel:
The integrator may be called through a .pbs file if operating on a scheduler.
