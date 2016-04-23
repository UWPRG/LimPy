# Langevin_Package (CURRENTLY OUTDATED README)
----
Langevin_Package is a Langevin Integrator capable of carrying out Molecular Dynamics
(MD), Metadynamics, Well-Tempered Metadynamics, and Infrequent Metadynamics.
The integrator is capable of operating on, and visualizing, 1-D and 2-D potentials.
This package can be operated to run from the command line (`run_langevin.py`) or through a GUI (`run_GUI.py`). A version of this code has also been created with MPI4Py so it may be implemented in parallel on multiple nodes (`run_langevin_parallel.py`)

----
## Use:

----
### Inputs:
**Input File**
The inputs specifying the system, method, parameters, and visualization options can be received via an inputs file or through a user interface.

The inputs files contains all of the specified inputs for the system of interest in the format below. Note information following #'s are comments indicating other options for those inputs.

```
Parameter,Input
Method,Well-Tempered Metadynamics
Potential_Function,cosine_potential
Plotting,False
Steps,100000
Step size,0.01
Temperature,300
Mass,1
X0,3.14
Xmin,-6.28
Xmax,6.28
Xincrement,0.01
Y0,0.5
Ymin,-2
Ymax,2
Yincrement,0.01
Gaussian Height,0.25
Gaussian Width,0.25
Deposition Frequency,500
Well Temperature,10000
Trials,100
Data Filename,PotclassTest_10000stride
Gamma,5
Kb,0.001987
Plot Freq,1000
Make Movie,False
X Boundary Condition,Quartic
Upper X BC,True
Lower X BC,False
Potential Parameters,2.5 0 2.5 0
```
**All lines must be included for a properly formatted input file, however a line is not applicable to a system (i.e. Y parameters for a 1-D system), then just ignore those contents or change the values to 0.**

Each of the definitions can be defined in the key below.

----
 ![alt text](https://github.com/UWPRG/Chris_Scripts/blob/master/Langevin_Package/Images/Input_Definitions.png)
----
**GUI**
The user can also input parameters directly through a user interface. For the preloaded potential functions, default variables are preloaded as inputs! Just enter your input choices and push the Lock in Values button. The inputs are the very same for the ones in the input files, but can be brought up by clicking the HELP button.

The user interface also provides the option to read the inputs from an input file. Just make sure the input file is in the inputfiles folder, type the name in the textbox, check the option Read File for Inputs, and click Lock in Values.

![alt text](https://github.com/UWPRG/Chris_Scripts/blob/master/Langevin_Package/Images/GUI_Window.png)
----
### Adding a New Potential

While the Langevin_Package comes with potential functions already pre-defined, a user can add new ones to the package in just a few easy steps!

Open up `potential_class.py` and add a new function following the template below. Potential Functions are created as a class with preset functions and attributes. You can just create a child class of PotentialFunction1D or PotentialFunction2D and overwrite the functions as necessary. The outline of a PotentialFunction1D is given below.
```python
class PotentialFunction1D:
    """Basic Outline of a 1-D Potential Function"""
    dimension='1-D Potential'
    def __init__(self):
        self.parameters=np.array([0,0])
        self.rare_event = np.array([float("inf")*-1,float("inf")])

    def get_potential(self, coords):
        """
        Equation that returns the potential energy at a CV
        To be specified by child class
        """
        pass

    def set_parameters(self, values):
        """Overwrites default parameters for potential function"""
        self.parameters = values

    def get_force(self, coords):
        """
        Equation that returns the force at a CV
        To be specified by child class
        """
        pass

    def set_lower_rare_event(self, lower):
        """Overwrites default lower rare event for potential function"""
        self.rare_event[0] = lower

    def set_upper_rare_event(self, upper):
        """Overwrites default upper rare event for potential function"""
        self.rare_event[1] = upper

    def get_triggered(self, coords):
        """
        Determines if rare event has happened, and specifies which one if
        there are multiple
        """
        trigger = False
        event = "null"
        if coords < self.rare_event[0]:
            trigger = True
            event = "A"
        elif coords > self.rare_event[1]:
            trigger = True
            event = "B"
        return trigger, event
```

----
## More to Come:
----
- **Selectivity**

- **Implement Reflective Boundary Conditions**

- **Unit Tests**

- **Bootstrapping of Error Analysis**
