# Python-Solver

A versatile python solver used for solving generic nonlinear systems. Developed as part of a collaboration, this is the solver part of it.

This is meant for solving systems of equations that can be written in the form R(U) = F(U,U)+S(U) = 0, sums of fluxes and sources. It uses a control volume like approach with more abstraction, relying on dictionaries to connect control volumes, allowing for more versatility, such as control volumes having different state variables, fluxes being one directional, and so forth. The underlying solver is scipy.optimize.fsolve with the remaining code handling construction and assembly of the system.



There are four files:
* [blocks.py] - contains definitions of block objects, our control volume like object
* [flux.py] - contains definitions of fluxes and flux functions
* [source.py] - contains definitions of sources and source functions
* [problem.py] - contains solvers, manages system construction

### Dependencies

The only optional dependency is numdifftools, used to obtain a numerical Jacobian of the system. If not already installed, it can be obtained with

pip install numdifftools

### Documentation 

Documentation is in doc/guide.tex. An example is in doc/examples.tex, with code in example.py. 

### Examples

Two other examples exist, using it to solve Poisson's equation and a diffusion equation in 2D, in poisson2D.py and diffusion2D.py, effectively using this as a PDE solver. In practice, it was developed for more arbitrary, smaller systems of equations.

### Tests

Tests are run with python test.py, which runs both poisson2D.py and diffusion2D.py and checks the result.