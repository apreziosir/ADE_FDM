#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D Transport equation solver for FDM with FTCS scheme.
High order derivatives are applied in the points where it is possible and an 
upwinding scheme is implemented for advection. 
The temporal term is treated with a second order Adams-Bashforth scheme
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana
February 2018
"""
# ==============================================================================
# Importing libraries and external functions
# ==============================================================================

import numpy as np  
import matplotlib.pyplot as plt
import time
import Screen_msgs as SM
import Analyt as AN
import Auxiliary as AUX

# ==============================================================================
# Declaration of physical variables for the problem. A positive value in X means
# going to the right and a positive value in Y means going up
# ==============================================================================

X0 = 0.                                         # Initial x point (m)
XL = 5.                                         # Final y point (m)
Y0 = 0.                                         # Initial y point (m)
YL = 5.                                         # Final y point (m)
Dx = 0.02                                       # Diff coeff x (m2/s)
Dy = 0.02                                       # Diff coeff y (m2/s)
u = 1.8                                         # Horizontal velocity (m/s)
v = 0.3                                         # Vertical velocity (m/s)
M = 10.                                         # Mass injected (g)
xC0 = 1.0                                       # x injection coordinate
yC0 = 1.0                                       # y injection coordinate
t0 = 1.                                         # Initial time (s)
A = (XL - X0) * (YL -Y0)                        # Domain area (m2)

# =============================================================================
# Numerical parameters for the program. For Crank-Nicholson ponderation, theta 
# factor means 0 = fully explicit, 1 = fully implicit
# =============================================================================

T = 10                                          # Total sim. time (s)
dt = 0.01                                      # timestep size (s)
Nx = 51                                        # Nodes in x direction
Ny = 51                                        # Nodes in y direction 
theta = 0                                       # Crank-Nicholson ponderation

# Calculation of initial parameters
nT = int(np.ceil((T - t0) / dt))                # Number of timesteps
dx = (XL - X0) / (Nx - 1)                       # dx (m)
dy = (YL - Y0) / (Ny - 1)                       # dy (m)

# Generating spatial mesh
x = np.linspace(X0, XL, Nx)                     # 1D array style 
y = np.linspace(Y0, YL, Ny)                     # 1D array style 
X, Y = np.meshgrid(x,y)                         # Meshgrid for plotting
del(x, y)                                       # Deleting useless arrays

# Generating error vectors (for svaing and comparing results)
errt = np.zeros(int(nT))                        # Error evolution in time

# ==============================================================================
# Boundary condtions for the program 1 = Dirichlet, 0 = Neumann.
# BC are stored in the following order [Bottom, Left, Right, Top] 
# TBC = Type of boundary condition 
# VBC = Value of boundary condition (g/m2) if Dirichlet (g/m2/m) if Neumann
# BBi, RBi, LBi, TBi identify de node numbers that are part of the domain
# ==============================================================================

TBC = np.array([0, 1, 0, 1])
VBC = np.array([0, 0, 0, 0])

BBi, LBi, RBi, TBi = AUX.boundaries(Nx, Ny)
ORi = AUX.ring(Nx, Ny)

# ==============================================================================
# Calculation of nondimensional parameters of the ADE (Sx, Sy, CFLx and CFLy)
# ==============================================================================

Sx = Dx * dt / (dx ** 2)
Sy = Dy * dt / (dy ** 2)
CFLx = u * dt / dx
CFLy = v * dt / dy

SM.show_sc(Sx, Sy, CFLx, CFLy)

# ==============================================================================
# Imposing initial condition as an early analytical solution of the ADE equation
# It is C(x0, y0, t0).
# ==============================================================================

C0 = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0)

# Estimating Cmax for plotting purposes
Cmax = np.max(C0)

# ==============================================================================
# First time step (must be done in a forward time scheme since there are no 
# other points in time). The first step has a size dt, hence the total time is 
# 1 * dt
# ==============================================================================

# Calculating analytical solution for the first timestep
Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + dt)

C1e = np.zeros((Nx, Ny))

# Imposing boundary conditions - matrix style (EXPLICIT METHOD)
C1e[0, :] = Ca[0, :] 
C1e[Nx - 1, :] = Ca[Nx -1, :]
C1e[:, 0] = Ca[:, 0]
C1e[:, Ny - 1] = Ca[:, Ny -1] 

# Solving the outer ring - matrix style (low order derivatives) (EXPLICIT 
# METHOD)

# Bottom part of the ring
for i in range(1, Nx - 1):
    
    C1e[1, i] = C0[1, i + 1] * Sx + C0[1, i - 1] * (Sx + CFLx) + C0[2, i] * \
    Sy + C0[0, i] * (Sy + CFLy) + C0[1, i] * (1 - 2 * Sx - 2 * Sy - CFLx - CFLy)

# Left part of the ring
for i in range(2, Ny - 1):
    
    C1e[i, 1] = C0[i + 1, 1] * Sy + C0[i - 1, 1] * (Sy + CFLy) + C0[i, 0] * \
    (Sy + CFLy) + C0[i, 2] * Sy + C0[i, 1] * (1 - 2 * Sx - 2 * Sy - CFLx - CFLy)

# Right part of the ring
for i in range(2, Ny - 1):
    C1e[i, Nx - 2] = C0[i + 1, Nx - 2] * Sy + C0[i - 1, Nx - 2] * (CFLy + Sy) \
    + C0[i, Nx - 1] * Sx + C0[i, Nx - 3] * (Sx + CFLx) + C0[i, Nx - 2] * \
    (1 - 2 * Sx - 2 * Sy - CFLx - CFLy)

# Top part of the ring
for i in range(2, Nx - 2):
    C1e[Ny - 2, i] = C0[Ny - 1, i] * Sy + C0[Ny - 3, i] * (Sy + CFLy) + \
    C0[Ny - 2, i + 1] * Sx + C0[Ny - 2, i - 1] * (Sx + CFLx) + C0[Ny - 2, i] * \
    (1 - 2 * Sx - 2 * Sy - CFLx - CFLy)

# Solving inner portion of the domain with forward time derivative and high 
# order spatial terms (just for first time step)
for i in range(2, Ny - 2):
    
    for j in range(2, Nx - 2):

#        # Low order derivatives (just for testing purposes, though it can be 
#        # implemented)
#        C1e[i, j] = Sx * C0[i, j + 1] + (Sx + CFLx) * C0[i, j - 1] + Sy * \
#        C0[i + 1, j] + (Sy + CFLy) * C0[i - 1, j] + C0[i, j] * (1 - 2 * Sx - \
#        2 * Sy - CFLx - CFLy)

        # High order approximation for implementation
        C1e[i, j] = -C0[i + 2, j] * (Sy / 12) + C0[i + 1, j] * (4 * Sy / 3) + \
        C0[i - 1, j] * (2 * CFLy + 4 * Sy / 3) - C0[i - 2, j] * (0.5 * CFLy + \
        Sy / 12) - C0[i, j + 2] * (Sx / 12) + C0[i, j + 1] * (4 * Sx / 3) + \
        C0[i, j - 1] * (2 * CFLx + 4 * Sx / 3) - C0[i, j - 2] *(Sx / 12 + \
        0.5 * CFLx) + C0[i , j] * (1 - 5 * Sx / 2 - 5 * Sy / 2 - 3 * CFLx / 2 \
        - 3 * CFLy / 2)
        
# Checking error and saving it in the assigned array
errt[0] = np.linalg.norm(C1e - Ca)

# ==============================================================================
# General time loop after first time step (implementing the leapfrog explicit
# scheme)
# ==============================================================================

# Defining vector for t+1 step
C2e = np.zeros((Nx, Ny))

# Activating interactive plotting and creating figure for interactive purpose of
# the program
plt.ion()
fig, axarr = plt.subplots(2, 2)

# Defining plots
cf1 = axarr[0, 0]
cf2 = axarr[0, 1]
cf3 = axarr[1, 0]
cf4 = axarr[1, 1]

fig.colorbar(cf1, ax=axarr[0, 0], ticks = np.linspace(0, Cmax, 5))

for I in range(2, 100):
    
    # Calculating analytical solution for the given timestep
    Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + I * dt)
    
    # Imposing boundary conditions
    C2e[0, :] = Ca[0, :] 
    C2e[Nx - 1, :] = Ca[Nx -1, :]
    C2e[:, 0] = Ca[:, 0]
    C2e[:, Ny - 1] = Ca[:, Ny -1] 
    
    # Calculating outer ring (low order derivatives)
    # Bottom part of the ring AB IMPLEMENTED
    for j in range(1, Nx - 1):
        
        C2e[1, j] = C1e[1, j + 1] * (2 * Sx / 3) + C1e[1, j - 1] * (2 / 3) * \
        (Sx + CFLx) + C1e[2, j] * (2 * Sy / 3) + C1e[0, j] * (2 / 3) * (Sy + \
        CFLy) + C1e[1, j] * (1 / 3) * (-4 * Sx - 4 * Sy - 2 * CFLx - 2 * CFLy) \
        - C0[1, j] * (1 / 3)

    # Left part of the ring AB IMPLEMENTED
    for j in range(2, Ny - 1):
        
        C2e[j, 1] = C1e[j + 1, 1] * (2 * Sy / 3) + C1e[j - 1, 1] * (2 / 3) * \
        (Sy + CFLy) + C1e[j, 0] * (2 / 3) * (Sy + CFLy) + C1e[j, 2] * (2 * Sy \
        / 3) + C1e[j, 1] * (1 / 3) * (-4 * Sx - 4 * Sy - 2 * CFLx - 2 * CFLy) \
        - C0[j, 1] * (1 / 3)
        
    # Right part of the ring AB IMPLEMENTED
    for j in range(2, Ny - 1):
        
        C2e[j, Nx - 2] = C1e[j + 1, Nx - 2] * (2 * Sy / 3) + C1e[j - 1, \
        Nx - 2] * (2 / 3) * (CFLy + Sy) + C1e[j, Nx - 1] * (2 * Sx / 3) + \
        C1e[j, Nx - 3] * (2 / 3) * (Sx + CFLx) + C1e[j, Nx - 2] * (1 / 3) * \
        (-4 * Sx - 4 * Sy - 2 * CFLx - 2 * CFLy) - C0[j, Nx - 2] * (1 / 3)

    # Top part of the ring AB IMPLEMENTED
    for j in range(2, Nx - 2):
        
        C2e[Ny - 2, j] = C1e[Ny - 1, j] * (2 * Sy / 3) + C1e[Ny - 3, j] * \
        (2 / 3) * (Sy + CFLy) + C1e[Ny - 2, j + 1] * (2 * Sx / 3) + \
        C1e[Ny - 2, j - 1] * (2 / 3) * (Sx + CFLx) + C1e[Ny - 2, j] * (1 / 3) \
        * (-4 * Sx - 4 * Sy - 2 * CFLx - 2 * CFLy) - C0[Ny - 2, j] * (1 / 3)
    
    # Calculating inner nodes (high order derivatives)
    for i in range(2, Ny - 2):
    
        for j in range(2, Nx - 2):

#            # Low order derivatives (just for testing purposes, though it can 
#            # be implemented)
#            C2e[i, j] = Sx * C1e[i, j + 1] + (Sx + CFLx) * C1e[i, j - 1] + Sy \
#            * C1e[i + 1, j] + (Sy + CFLy) * C1e[i - 1, j] + C1e[i, j] * (-2 * \
#            Sx - 2 * Sy - CFLx - CFLy) + C0[i, j]

            # High order approximation for implementation
            C2e[i, j] = -C1e[i + 2, j] * (Sy / 18) + C1e[i + 1, j] * (8 * Sy / \
            9) + C1e[i - 1, j] * (4 * CFLy / 3 + 8 * Sy / 9) - C1e[i - 2, j] * \
            (CFLy / 3 + Sy / 18) - C1e[i, j + 2] * (Sx / 6) + C1e[i, j + 1] * \
            (8 * Sx / 9) + C1e[i, j - 1] * (4 * CFLx / 3 + 8 * Sx / 9) - \
            C1e[i, j - 2] * (Sx / 18 + CFLx / 3) + C1e[i , j] * (-5 * Sx / 3 - \
            5 * Sy / 3 - CFLx - CFLy + 4 / 3) - C0[i, j] * (1 / 3)
    
    
    # Storing error in the errt vector
    errt[I - 1] = np.linalg.norm(C2e - Ca)
    
    # Plotting relevant data for the process
    cf = axarr[0 ,0].contourf(X,Y,C2e)
    axarr[0 ,0].set_title('Numerical solution')
    axarr[0, 1].contourf(X, Y, Ca)
    axarr[0, 1].set_title('Analytical solution')
    axarr[1, 0].contourf(X, Y, np.abs(C2e - Ca))
    axarr[1, 0].set_title('Error')
    plt.draw()
    plt.pause(0.1)
    
    # Changing the variables for the next step
    C0 = C1e
    C1e = C2e
    C2e = np.zeros((Nx, Ny))    