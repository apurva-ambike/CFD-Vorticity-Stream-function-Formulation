# CFD: Vorticity–Stream Function Formulation
Numerical solution of 2D incompressible flow using the vorticity–stream function formulation of the Navier–Stokes equations, applied to two problems: irrotational (potential) flow and viscous laminar flow over two surface-mounted square ribs.

## Overview

The vorticity–stream function approach recasts the incompressible Navier–Stokes equations into two coupled scalar PDEs — a Poisson equation for the stream function and a transport equation for vorticity — eliminating the pressure term and automatically satisfying continuity. This repository implements finite-difference solvers for two cases:

- **Irrotational (potential) flow** — an inviscid benchmark consisting of a stream function Laplace equation used to validate the linear solvers and grid setup.
- **Viscous flow over surface-mounted ribs** — 2D laminar flow over two surface-mounted square ribs, including boundary layer development, recirculation zones and wake formation.

Several iterative solvers to evaluate systems of linear equations are implemented to solve the equations associated with the two problem cases, including Gauss–Seidel, Successive Over-Relaxation (SOR), BiCGSTAB, BiCGSTAB with ILU(0) preconditioning, and BiCGSTAB with SIP (Strongly Implicit Procedure) preconditioning. Their convergence behaviour is analysed and the effect of preconditioning is studied. 

## Repository Structure

```
├── Grid_generator/          # Mesh/grid generation codes
├── Irrotational Case/       # Inviscid potential flow solved using Gauss-Seidel, SOR, BiCGSTAB and BiCGSTAB with ILU(0) and SIP preconditioning
├── Viscous case/            # Laminar viscous flow solver (vorticity transport equation and stream function Poisson's equation, each solved using BICGSTAB + SIP)
├── Plotting_python_files/   # Post-processing and visualization scripts (matplotlib)
├── Report & Documents/      # The accompanying technical report
├── Report_images/           # Figures used in the report
├── Results/                 # Simulation output data in .csv format
└── README.md
```
## Requirements

- C++ compiler (C++11 or later recommended)
- Python 3.14 or later versions with `numpy` and `matplotlib` for post-processing

## Notes

1) The .csv file which contains the x-coordinates for the grid is named xcoords_3.csv on purpose (there are no files like xcoords_1.csv or xcoords_2.csv) as the domain length was extended from its intial span. Similarly, the .txt file required for generating the x-coordinates is named grid_2.txt on purpose.
2) The code for the BiCGSTAB solver without preconditioning (BiCGSTAB.cpp) contains functions returning and using long double data type instead of the regular double data type. This is because using double data type, results in it being unable to converge to the decided tolerance of 1e-8 (rather, it starts to stagnate around 1.5e-8)
   
## Acknowledgements

This work was carried out during a research internship at IIT Guwahati under the guidance of Prof. Arnab Kr. De.

## Author

Apurva Ambike
