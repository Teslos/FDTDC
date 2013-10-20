# Usage
The program quantfdtd implements FDTD method to solve Schrodinger equation of the electron
in arbitary potential.

#Input file
Input file defines the variables and parameters used in the calculations.
```
# Number of cells in x direction
CELLSX=100
# Number of cells in y directions
CELLSY=100
# Size of cell in x direction
SIZEX=1.e-10
# conversion factor
CONV=1e9
# timestep
TIMESTEP=0.04e-16
# electronic mass
EMASS=0.42
# potential
POTENTIAL=6
# number of steps
NSTEPS=1200000
# eigenvalue of quantum dot
EIGEN=nan
# define the test point 
TESTX=50
TESTY=50
```
 
