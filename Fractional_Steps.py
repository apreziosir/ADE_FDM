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
import scipy.sparse as scsp
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
v = 0.9                                        # Vertical velocity (m/s)
M = 1.                                          # Mass injected (g)
xC0 = 0.0                                       # x injection coordinate
yC0 = 0.0                                       # y injection coordinate
t0 = 1.                                         # Initial time (s)
A = (XL - X0) * (YL -Y0)                        # Domain area (m2)

# =============================================================================
# Numerical parameters for the program. 
# =============================================================================

T = 3                                           # Total sim. time (s)
dt = 0.01                                       # timestep size (s)
Nx = 21                                         # Nodes in x direction
Ny = 21                                         # Nodes in y direction 

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



