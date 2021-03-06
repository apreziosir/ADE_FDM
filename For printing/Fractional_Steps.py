#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scritp dor calculating ADE equation with fractional steps
Advecyion is solved with Adams Bashforth second order approach and the diffusive 
term with a Crank Nicholson scheme
Created on Mon Mar  5 15:38:07 2018
@author: Antonio Preziosi-Ribero
CFD Pontificia Universidad Javeriana
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib import cm
from scipy.sparse.linalg import spsolve
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
Dx = 0.030                                      # Diff coeff x (m2/s)
Dy = 0.045                                      # Diff coeff y (m2/s)
u = 0.8                                         # Horizontal velocity (m/s)
v = 0.9                                         # Vertical velocity (m/s)
M = 1.                                          # Mass injected (g)
xC0 = 0.0                                       # x injection coordinate
yC0 = 0.0                                       # y injection coordinate
t0 = 1.                                         # Initial time (s)
A = (XL - X0) * (YL -Y0)                        # Domain area (m2)

# =============================================================================
# Numerical parameters for the program. 
# theta = 1 fully implicit
# theta = 0 fully explicit
# =============================================================================

T = 3                                           # Total sim. time (s)
dt = 0.001                                       # timestep size (s)
Nx = 21                                         # Nodes in x direction
Ny = 21                                         # Nodes in y direction 
theta = 0.5                                     # Crank-Nicholson ponderation

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
# Calculation of nondimensional parameters of the ADE (Sx, Sy, CFLx and CFLy)
# ==============================================================================

Sx = Dx * dt / (dx ** 2)
Sy = Dy * dt / (dy ** 2)
CFLx = u * dt / dx
CFLy = v * dt / dy

SM.show_sc(Sx, Sy, CFLx, CFLy)

# Matrix assembly for diffusive implicit part
LHS = AUX.Mat_ass(Sx, Sy, Nx, Ny)

# ==============================================================================
# Imposing initial condition as an early analytical solution of the ADE equation
# It is C(x0, y0, t0).
# ==============================================================================

C0 = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0)

# Estimating Cmax for plotting purposes
Cmax = np.max(C0)

# ==============================================================================
# First time step. Forward Euler for advective part, and Crank-Nicholson for 
# diffusive term
# ==============================================================================

# First fractional step part for timestep 1 - advective term

# Calculating analytical solution for first timestep
Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + dt)

# Evaluating space derivative just for advections
sp0 = AUX.ev_sp(C0, 0, 0, dx, dy, u, v, Nx, Ny)

# Making first advection step
C1s = C0 + dt * sp0

# Imposing boundary conditions (domain is a 2D array)

C1s[0, :] = Ca[0, :]                        # Bottom boundary
C1s[Ny - 1, :] = Ca[Ny -1, :]               # Top boundary
C1s[:, 0] = Ca[:, 0]                        # Left boundary
C1s[:, Nx - 1] = Ca[:, Nx - 1]              # Right boundary

# Second fractional step part for timestep 1 - diffusive term

# Explicit part - evaluating diffusion just in space
sp0 = AUX.ev_sp(C1s, Dx, Dy, dx, dy, 0, 0, Nx, Ny)

C1e = C1s + dt * sp0

# Implicit part
C1i = spsolve(LHS, np.reshape(C1s, (Nx * Ny, 1)))

C1 = theta * C1i.reshape(Nx, Ny) + (1 - theta) * C1e

# Error checking
errt[0] = np.linalg.norm(C1 - Ca)

# Second array declaration
C2s = np.zeros((Ny, Nx))
err = np.zeros((Ny, Nx))

# Figure declaration 
plt.ion()
plt.figure(1, figsize=(14,8))
#style.use('fivethirtyeight')
style.use('ggplot')


# Starting time loop
for I in range(2, nT + 1):
    
    # Calculating analytical solution for the time step
    Ca = AN.difuana(M, A, Dx, Dy, u, v, xC0, yC0, X, Y, t0 + I * dt)
    
    # Calculating C* - only advection with Adams Bashforth Order 2
    
    # Spatial function evaluated in t - 1
    sp0 = AUX.ev_sp(C0, 0, 0, dx, dy, u, v, Nx, Ny)
    
    # Spatial function evaluated in t 
    sp1 = AUX.ev_sp(C1, 0, 0, dx, dy, u, v, Nx, Ny)
    
    # Implementing 2nd order Adams-Bashforth
    C2s = C1 + (dt / 2) * (3 * sp0 - sp1)
    
    # Imposing Boundary conditions
    C2s[0, :] = Ca[0, :] 
    C2s[Ny - 1, :] = Ca[Ny -1, :]
    C2s[:, 0] = Ca[:, 0]
    C2s[:, Nx - 1] = Ca[:, Nx - 1]
    
    # Calculating C2 with the values of C* - explicit FE
    C2sp = AUX.ev_sp(C2s, Dx, Dy, dx, dy, 0, 0, Nx, Ny)
    C2e = C2s + dt * C2sp
    
    # Calculating C2 with C* - implicit BE
    C2i = spsolve(LHS, np.reshape(C2s, (Nx * Ny, 1)))
    
    # Making Crank-Nicholson ponderation
    C2 = theta * C2i.reshape(Nx, Ny) + (1 - theta) * C2e
    
    # Calculating error for the timestep given
    err = np.abs(C2 - Ca)
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
    C0 = C1
    C1 = C2
    C2 = np.zeros((Ny, Nx))