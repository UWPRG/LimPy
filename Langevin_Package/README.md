# Langevin_Package (Still a work in progress)
----
Langevin_Package is a Langevin Integrator capable of carrying out Molecular Dynamics
(MD), Metadynamics, Well-Tempered Metadynamics, and Infrequent Metadynamics.
The integrator is capable of operating on, and visualizing, 1-D and 2-D potentials.
This package can be operated to run from the command line (`run_langevin.py`) or through a GUI (`run_langevin_GUI.py`). A version of this code has also been created with MPI4Py so it may be implemented in parallel on multiple nodes (`run_langevin_parallel.py`)

----
## Use:
----
### Inputs:
**Input File**
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
**All lines must be included for a properly formatted input file, however a line is not applicable to a system (i.e. Y parameters for a 1-D system), then just ignore those contents or change the values to 0.**

Each of the definitions can be defined in the key below.

----
 ![alt text](https://github.com/UWPRG/Chris_Scripts/blob/master/Langevin_Package/Images/Input%20Definitions.png)
----
**GUI**
The user can also input parameters directly through a user interface. For the preloaded potential functions, default variables are preloaded as inputs! Just enter your input choices and push the Lock in Values button. The inputs are the very same for the ones in the input files, but can be brought up by clicking the HELP button.

The user interface also provides the option to read the inputs from an input file. Just make sure the input file is in the inputfiles folder, type the name in the textbox, check the option Read File for Inputs, and click Lock in Values.

![alt text](https://github.com/UWPRG/Chris_Scripts/blob/master/Langevin_Package/Images/GUI_window.png)
----
### Adding a New Potential

While the Langevin_Package comes with potential functions already pre-defined, a user can add new ones to the package in just a few easy steps!

1) Open up `potential_functions.py` and add a new function following the template below. Just fill in the inputs.
```python
def your_potential(coords):
    """
    Calculate the Potential Energy and Force. Determine if rare event has happened.

    Parameters:
    -----------

    coords  : float (if 1-D) or array of floats
              If 2-D the replace coords with (x,y)


    Returns:
    --------

    V       : float (or array of floats)
              Potential Energy

    F       : float (or array of floats)
              Force

    Trigger : Boolean
              Has rare event occurred (True) or not (False)
    """
    V = # Input Potential Equation
    F = # Input Force Equation

    if (# Input Rare Event
        and if hasattr(coords, "__len__") is True):
        Trigger = True
    else:
        Trigger = False

    return (V, F, Trigger)
```
2) Then add the new function to the dictionary of function in `get_potential_dict()` in `potential_functions.py`.

3) Define the boundary condition function for the new function (even if there is none).

```python
def your_potential_bc(vnew,f2,coords):
    """
    Apply Boundary Condition to the potential, force, and coordinates.

    Parameters:
    -----------
    vnew       : float (or array of floats)
                 Potential Energy

    f2         : float (or array of floats)
                 Force

    coords     : float
                 coordinates

    Returns:
    --------

    vnew       : float (or array of floats)
                 Adjusted potential energy from boundary condition

    F          : float (or array of floats)
                 Adjusted force from boundary condition

    coords     : float
                 adjusted coordinates from boundary condition

    bcbias     : float
                 bias applied strictly from the boundary condition

    is_periodic: Boolean
                 indicates if the boundary is periodic (True) or not (False)
    """
    # Add changes that need to occur to potential, force, coordinates, or bias.
    return (vnew, f2, coords, bcbias, is_periodic)
```
4) Add the new function to `get_boundary_condition_dict()` in `potential_functions.py` with the same key as the entry in the `get_potential_dict()`

5) You can now call your new function from an input file using the key for the dictionary entries or select it from the dropdown menu in the GUI.
