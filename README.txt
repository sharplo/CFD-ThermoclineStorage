This program consists of 10 files: Input.f95, Measure.f95, Dynamics.f95, OrderVerif.f95, Cycle.f95, main.f95, plot.py, makefile, parameters.dat, sol-exact.dat

InOut.f95 declares global variables and constants, and contains subroutines that write data onto files
Measure.f95 contains functions that returns measured quantities such as error norms
Dynamics.f95 contains subroutines that evolve the approximate solutions
OrderVerif.f95 contains subroutines that run procedures of order verification study
Cycle.f95 contains subroutines that run simulation in storage cycles
main.f95 reads in parameters from NAMELIST and governs which of the studies to be carried out
plot.py plots the data obtained from the simultaions and saves them in a folder called "figures"
makefile compiles the fortran source codes and executes the programs
parameters.dat contains NAMELIST parameters that the program use
sol-exact.dat contains exact solutions for comparison

To compile the fortran source codes, run the command "make" in which gfortran compiler is required.
Before executing the program, parameters should be specified in "parameters.dat".
To execute the program, either run the executable "thermocline" or run the command "make run" which will run the simulation and plot the graphs in which gfortran and python3 compiler are required. When the program is executed in the first time, the command "export OMP_NUM_THREADS = <number>" is required if not specified before. Also, the procedures for plotting in the bottom of plot.py can be commented out becasue the .dat files do not exist.

results.ods is a spreadsheet that contains parameters for different studies that can reproduce the results in the project.

Compiler used: GNU Fortran 7.4.0 and Python 3.6.9

Lo Chim Yui
10 December, 2019
