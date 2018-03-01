#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auxiliary functions for the 2D transport equation 
@author: Antonio Preziosi-Ribero
CFD- Pontificia Universidad Javeriana
"""

import numpy as np

# ==============================================================================
# Function tha identify the nodes that are boundary for the 2D transport 
# BBi = Bottom, LBi = left boundary, RBi = right boundary, TBi = top boundary
# ==============================================================================

def boundaries(Nx, Ny):
    
    BBi = np.linspace(0, Nx - 1, Nx)
    LBi = np.linspace(Nx, Nx * (Ny - 2), Ny - 2)
    RBi = np.linspace(2 * Nx - 1, Nx * (Ny - 1) - 1, Ny - 2)
    TBi = np.linspace(Nx * (Ny - 1), Nx * Ny - 1, Nx)
    
    return BBi, LBi, RBi, TBi

# ==============================================================================
# Function that identifies the outer ring where the high order derivatives 
# cannot be implemented (OR stands for outer ring) 
# ==============================================================================
    
def ring(Nx, Ny):
    
    BR = np.linspace(Nx + 1, 2 * Nx - 2, Nx -2)
    LR = np.linspace(Nx * 2 + 1, Nx * (Ny - 3) + 1, Ny - 4)
    RR = np.linspace(Nx * 3 - 2, Nx * (Ny - 2) - 2, Ny - 4)
    TR = np.linspace(Nx * (Ny -2) + 1, Nx * (Ny - 1) - 2, Nx - 2)
    OR = np.sort(np.concatenate((BR, LR, RR, TR), axis = 0))
    
    return OR

# ==============================================================================
# Evaluation space functions in outer ring for the Adams Bashforth scheme
# Inputs are the concentration vector and the numerical parameters Sx, Sy, CFLx, 
# and CFLy
# ==============================================================================
    
def ev_sp(C, Dx, Dy, dx, dy, u, v, Nx, Ny):
    
    # Boundaries are not evaluated in this function since they are imposed in 
    # the main routine of the program. Also, this function works for both the 
    # first time step and the subsequent steps since it just analyzes the value
    # of the spatial functions applied to C. 
    
    sp = np.zeros((Ny, Nx))
    
    # Bottom part of the ring
    for i in range(1, Nx - 1):
        
        sp[1, i] = C[2, i] * (Dy / dy ** 2) + C[0, i] * (v / dy + Dy / dy ** \
        2) + C[1, i + 1] * (Dx / dx ** 2) + C[1, i - 1] * (u / dx + Dx / dx ** \
        2) - C[1, i] * (u / dx + v / dy + 2 * Dx / dx ** 2 + 2 * Dy / dy ** 2)
    
    # Left part of the ring    
    for i in range(1, Ny - 1):
        sp[i, 1] = C[i + 1, 1] * (Dy / dy ** 2) + C[i - 1, 1] * (v / dy + Dy \
        / dy ** 2) + C[i, 2] * (Dx / dx ** 2) + C[i, 0] * (u / dx + Dx / dx ** \
        2) - C[i, 1] * (u / dx + v / dy + 2 * Dx / dx ** 2 + 2 * Dy / dy ** 2)
    
    # Right part of the ring    
    for i in range(2, Ny - 1):
        sp[i, Nx - 2] = C[i + 1, Nx - 2] * (Dy / dy ** 2) + C[i - 1, Nx - 2] \
        * (v / dy + Dy / dy ** 2) + C[i, Nx - 1] * (Dx / dx ** 2) + C[i, Nx - \
        3] * (u / dx + Dx / dx ** 2) - C[i, Nx - 2] * (u / dx + v / dy + 2 * \
        Dx / dx ** 2 + 2 * Dy / dy ** 2)
            
    # Top part of the ring
    for i in range(2, Nx - 2):
        
        sp[Ny - 2, i] = C[Ny - 1, i] * (Dy / dy ** 2) + C[Ny - 3, i] * (v / dy \
        + Dy / dy ** 2) + C[Ny - 2, i + 1] * (Dx / dx ** 2) + C[Ny - 2, i - 1] \
        * (u / dx + Dx / dx ** 2) - C[Ny - 2, i] * (u / dx + v / dy + 2 * Dx / \
        dx ** 2 + 2 * Dy / dy ** 2)
        
    # Inner portion of the domain
    for i in range(2, Ny - 2):
#    
        for j in range(2, Nx - 2):
            
            sp[i ,j] = C[i + 1, j] * (Dy / dy ** 2) + C[i - 1, j] * (Dy / dy \
            ** 2 + 2 * v / dy) - C[i - 2, j] * (v / (2 * dy)) + C[i, j + 1] * \
            (Dx / dx ** 2) + C[i, j - 1] * (Dx / dx ** 2 + 2 * u / dx) - C[i, \
            j - 2] * (u / (2 * dx)) - C[i, j] * (2 * Dx / dx ** 2 + 2 * Dy / \
            dy ** 2 + (3 * v) / (2 * dy) + (3 * u) / (2 * dx))
        
    return sp
        
    