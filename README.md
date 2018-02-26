# ADE_FDM

Solution to the two dimensional advection dispersion equation with high order finite differences with explicit-implicit scheme. 
The time discretization is made by implementing a Leapfrog scheme, high order upwinding for the advective term and high order second derivative for the dispersion term. This scheme is compared with the analytical solution of the equation given in the Analytical.py script.

Boundary conditions can be modified to whatever value of Neumann or Dirichlet. However, for the sake of the example, the BC imposed in all the domain boundary is the analytical solution presented in Analytical.py

Afterwards, the idea is to implement a fractional step method and compare its results to the high order finite differences method. 
