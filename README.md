# NeutronSpinTransport
This is a simulation written in C++ to track magnetization of individual neutrons in a beam through a changing magnetic field.


# Program Structure
Neutron Generator (generate.hpp and generate.cpp) generates neutrons within a 6cm circular aperture with a speed distribution identical to a 
Maxwellian distribution. Generated neutrons are then followed down the apparatus using the values of angular deviations and velocities. 

Runge-Kutta Solver (RK4Solver.h, RK4Solver.cpp, DESolver.h, DESolver,cpp, DiffEqn.h and DiffEqn.h) was used to solve a set of differential equations (Bloch equation) in order to calculate the 
values of magnetization of individual neutrons.

Magnetic Field (Bfield.hpp and Bfield.cpp) used in the simulation was provided to the solver in order to perform the numerical analysis.

Polarization (track.h and track.cpp) of the neutron beam was found using the RK-Solver and the given magnetic field.

This simulation outputs the data in csv files which are then fed to MATLAB to do further visual and quantitative analysis.
