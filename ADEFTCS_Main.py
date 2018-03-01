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
import matplotlib.animation as animation
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
u = 0.3                                         # Horizontal velocity (m/s)
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
dt = 0.1                                        # timestep size (s)
Nx = 11                                         # Nodes in x direction
Ny = 11                                         # Nodes in y direction 

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

C1e = np.zeros((Ny, Nx))

C1e = C0 + dt * AUX.ev_sp(C0, Dx, Dy, dx, dy, u, v, Nx, Ny)

# Imposing boundary conditions - matrix style (EXPLICIT METHOD)
C1e[0, :] = Ca[0, :] 
C1e[Ny - 1, :] = Ca[Ny -1, :]
C1e[:, 0] = Ca[:, 0]
C1e[:, Nx - 1] = Ca[:, Nx - 1] 

errt[0] = np.linalg.norm(C1e - Ca)

# Entering the main time loop and calculating each concentration field 
# numerically and analytically

# Second array declaration
C2e = np.zeros((Ny, Nx))
err = np.zeros((Ny, Nx))

# Initiating plots
plt.ion()
style.use('fivethirtyeight')
fig = plt.figure()
fig.set_size_inches(14, 12)
fig.subplots_adjust(hspace=0.15)
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)


# ==============================================================================
# Defining plotting for animation function
# ==============================================================================

#def animate(I):
#    
#    ax1.clear()
#    ax2.clear()
#    ax3.clear()
#    ax4.clear()
#    ax1.contourf(X, Y, C2e, cmap=cm.coolwarm)
#    ax1.set_title('Numerical solution')
#    ax2.contourf(X, Y, Ca, cmap=cm.coolwarm)
#    ax2.set_title('Analytical solution')
#    ax3.contourf(X, Y, err, cmap=cm.coolwarm)
#    ax3.set_title('Error field')   
#    ax4.semilogy(np.linspace(0, nT - 1, nT), errt)
#    ax4.set_title('Error in every timestep')
    
#Plotting part
#fig, axarr = plt.subplots(2, 2)
#fig.subplots_adjust(hspace=0.5)
#axarr[0, 0].contourf([], [], [], cmap = cm.coolwarm)
#axarr[0, 0].set_title('Numerical solution')
#fig.colorbar(axarr[0, 0], ax=axarr[0, 0], shrink=0.9)
#axarr[0, 1].contourf([], [], [], cmap = cm.coolwarm)
#axarr[0, 1].set_title('Analytical solution')
#axarr[1, 0].contourf([], [], [], cmap = cm.coolwarm)
#axarr[1, 0].set_title('Error heatmap')
#
#
#
#line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#                                   interval=50, blit=True)

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
    
    # Classical plotting - not working at all
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax1.contourf(X, Y, C2e, cmap=cm.coolwarm)
    ax1.set_title('Numerical solution')
    ax2.contourf(X, Y, Ca, cmap=cm.coolwarm)
    ax2.set_title('Analytical solution')
    ax3.contourf(X, Y, err, cmap=cm.coolwarm)
    ax3.set_title('Error field')   
    ax4.semilogy(np.linspace(0, nT - 1, nT), errt)
    ax4.set_title('Error in every timestep')
    fig.canvas.draw()
    
#    # Plotting the results in a window - working just for final step
#    ani = animation.FuncAnimation(fig, animate(I), interval = 1000, blit = True)
#    plt.show()
    
    # Updating concentration fields
    C0 = C1e
    C1e = C2e
    C2e = np.zeros((Ny, Nx))