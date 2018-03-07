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
from matplotlib import style
from matplotlib import cm
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
Dx = 0.30                                      # Diff coeff x (m2/s)
Dy = 0.85                                      # Diff coeff y (m2/s)
u = 1.2                                         # Horizontal velocity (m/s)
v = 0.9                                        # Vertical velocity (m/s)
M = 1.                                          # Mass injected (g)
xC0 = 0.0                                       # x injection coordinate
yC0 = 0.0                                       # y injection coordinate
t0 = 1.                                         # Initial time (s)
A = (XL - X0) * (YL -Y0)                        # Domain area (m2)

# =============================================================================
# Numerical parameters for the program. 
# =============================================================================

T = 4.                                           # Total sim. time (s)
dt = 0.0001                                       # timestep size (s)
Nx = 101                                         # Nodes in x direction
Ny = 101                                         # Nodes in y direction 

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
errt = np.zeros(int(nT))                       # Error evolution in time

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
# First time step (must be done in a forward Euler time scheme since there are 
# no  other points in time). The first step has a size dt, hence the total time 
# is (1 * dt)
# ==============================================================================

# Calculating analytical solution for the first timestep
Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + dt)

C1e = np.zeros((Ny, Nx))
sp0 = AUX.ev_sp(C0, Dx, Dy, dx, dy, u, v, Nx, Ny)

C1e = dt * sp0 + C0

# Imposing boundary conditions (domain is a 2D array)

C1e[0, :] = Ca[0, :]                        # Bottom boundary
C1e[Ny - 1, :] = Ca[Ny -1, :]               # Top boundary
C1e[:, 0] = Ca[:, 0]                        # Left boundary
C1e[:, Nx - 1] = Ca[:, Nx - 1]              # Right boundary


# Calculating first time step error
errt[0] = np.linalg.norm(C1e - Ca)

# Second array declaration
C2e = np.zeros((Ny, Nx))
err = np.zeros((Ny, Nx))

# Figure declaration 
plt.ion()
plt.figure(1, figsize=(14,8))
#style.use('fivethirtyeight')
style.use('ggplot')

# starting time loop
for I in range(2, nT + 1):
    
    # Calculating analytical solution for the time step
    Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + I * dt)
    
    # Spatial function evaluated in t
    sp0 = AUX.ev_sp(C0, Dx, Dy, dx, dy, u, v, Nx, Ny)
    
    # Spatial function evaluated in t - 1
    sp1 = AUX.ev_sp(C1e, Dx, Dy, dx, dy, u, v, Nx, Ny)
    
    # Implementing 2nd order Adams-Bashforth
    C2e = C1e + (dt / 2) * (3 * sp0 - sp1)
    
    # Imposing Boundary conditions
    C2e[0, :] = Ca[0, :] 
    C2e[Ny - 1, :] = Ca[Ny -1, :]
    C2e[:, 0] = Ca[:, 0]
    C2e[:, Nx - 1] = Ca[:, Nx - 1]
    
    # Calculating error for the timestep given
    err = np.abs(C2e - Ca)
    errt[I - 1] = np.linalg.norm(C2e - Ca)
    
    # Plotting numerical, analytical, error field and error in time
    plt.clf()
    
    # Numerical solution
    plt.subplot(2, 2, 1)
    plt.contourf(X, Y, C2e, cmap=cm.coolwarm, vmin = 0., vmax = Cmax)  
    plt.title('Numerical solution')     
    plt.colorbar()
    
    # Analytical solution
    plt.subplot(2, 2, 2)
    plt.contourf(X, Y, Ca, cmap=cm.coolwarm, vmin = 0., vmax = Cmax)                          
    plt.title('Analytical solution')
    plt.colorbar()
    
    # Error field
    plt.subplot(2, 2, 3)
    plt.contourf(X, Y, np.abs(err), cmap=cm.coolwarm)                         
    plt.title('Error field')
    plt.colorbar()
    
    # Error evolution
    plt.subplot(2, 2, 4)
    plt.semilogy(np.linspace(0, dt * (nT - 1), nT), errt)
    plt.title('Norm 2 Error evolution')
    plt.ylim([1e-8, 1e2])
    
    plt.draw()
    plt.pause(0.01)
    
    # Saving a sample for the work
    if dt * I == 2.:
        plt.savefig('Results_2s.pdf')
    
    # Updating concentration fields
    C0 = C1e
    C1e = C2e
    C2e = np.zeros((Ny, Nx))