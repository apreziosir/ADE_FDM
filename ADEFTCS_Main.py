#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D Transport equation solver for FDM with FTCS scheme.
High order derivatives are applied in the points where it is possible and an 
upwinding scheme is implemented for advection. 
The temporal term is treated with a Leap Frog approximation
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana
February 2018
"""
# ==============================================================================
# Importing libraries and external functions
# ==============================================================================

import numpy as np  
import matplotlib.pyplot as plt
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
Dx = 0.2                                        # Diff coeff x (m2/s)
Dy = 0.2                                        # Diff coeff y (m2/s)
u = 0.1                                         # Horizontal velocity (m/s)
v = 0.1                                         # Vertical velocity (m/s)
M = 4.                                          # Mass injected (g)
xC0 = 2                                         # x injection coordinate
yC0 = 2                                         # y injection coordinate
A = (XL - X0) * (YL -Y0)                        # Domain area (m2)

# =============================================================================
# Numerical parameters for the program. For Crank-Nicholson ponderation, theta 
# factor means 0 = fully explicit, 1 = fully implicit
# =============================================================================

T = 10                                          # Total sim. time (s)
dt = 0.1                                        # timestep size (s)
Nx = 11                                         # Nodes in x direction
Ny = 11                                         # Nodes in y direction 
theta = 0                                       # Crank-Nicholson ponderation

# Calculation of initial parameters
nT = np.ceil(T / dt)                            # Number of timesteps
dx = (XL - X0) / (Nx - 1)                       # dx (m)
dy = (YL - Y0) / (Ny - 1)                       # dy (m)

# Generating spatial mesh
x = np.linspace(X0, XL, Nx)                     # 1D array style 
y = np.linspace(Y0, YL, Ny)                     # 1D array style 
X, Y = np.meshgrid(x,y)                         # Meshgrid for plotting
del(x, y)                                       # Deleting useless arrays

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
# It is almost as if a pulse is injected in (x0, y0). The analytical solution 
# does not work with t = 0
# ==============================================================================

C0 = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, 1e-8)

# ==============================================================================
# First time step (must be done in a forward time scheme since there are no 
# other points in time). The first step has a size dt, hence the total time is 
# 1 * dt.
# ==============================================================================

# Calculating analytical solution for the first timestep
Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, dt)

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
    
    C1e[1, i] = C0[1, i + 1] * Sx + C0[1, i - 1] * (Sx + CFLx) + C0[2, i] * Sy + C0[0, i] * (Sy + CFLy) + C0[1, i] * (1 - 2 * Sx - 2 * Sy - CFLx - CFLy)
#
## Left part of the ring
#for i in range(2, Ny - 1):
#    C1e[i, 1] = 
#
## Right part of the ring
#for i in range(2: Ny - 1):
#    C1e[i, Nx - 2] =
#
## Top part of the ring
#for i in range(2:Nx - 2):
#    C1e[Ny - 2, i] =
#

# Solving the inner portion (high order aproximation) - matrix style (EXPLICIT 
# METHOD)