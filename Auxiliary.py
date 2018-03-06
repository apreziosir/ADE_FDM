#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auxiliary functions for the 2D transport equation 
@author: Antonio Preziosi-Ribero
CFD- Pontificia Universidad Javeriana
"""

import numpy as np
import scipy.sparse as scsp

# ==============================================================================
# Evaluation space functions for the explicit parts of proposed schemes
# Inputs are the concentration vector and the numerical parameters Dx, Dy, u, v, 
# dx and dy
# For the fractional time steps technique the coefficients that are not playing
# in the solution are equal to zero 
# ==============================================================================
    
def ev_sp(C, Dx, Dy, dx, dy, u, v, Nx, Ny):
    
    # Boundaries are not evaluated in this function since they are imposed in 
    # the main routine of the program. Also, this function works for both the 
    # first time step and the subsequent steps since it just analyzes the value
    # of the spatial functions applied to C. 
    
    sp = np.zeros((Ny, Nx))
    
    # Defining coefficients for filling spatial formula
    CDx = Dx / (dx ** 2)
    CDy = Dy / (dy ** 2)
    CAx = u / dx
    CAy = v / dy
    
    # Bottom part of the ring
    for i in range(1, Nx - 1):        
        sp[1, i] = C[2, i] * CDy + C[0, i] * (CAy + CDy) + C[1, i + 1] * (CDx) \
        + C[1, i - 1] * (CAx + CDx) - C[1, i] * (CAx + CAy + 2 * CDx + 2 * CDy)
    
    # Left part of the ring    
    for i in range(1, Ny - 1):
        sp[i, 1] = C[i + 1, 1] * CDy + C[i - 1, 1] * (CAy + CDy) + C[i, 2] * \
        CDx + C[i, 0] * (CAx + CDx) - C[i, 1] * (CAx + CAy + 2 * CDx + 2 * CDy)
    
    # Right part of the ring    
    for i in range(2, Ny - 1):
        sp[i, Nx - 2] = C[i + 1, Nx - 2] * CDy + C[i - 1, Nx - 2] * (CAy + \
        CDy) + C[i, Nx - 1] * CDx + C[i, Nx - 3] * (CAx + CDx) - C[i, Nx - 2] \
        * (CAx + CAy + 2 * CDx + 2 * CDy)
            
    # Top part of the ring
    for i in range(2, Nx - 2):        
        sp[Ny - 2, i] = C[Ny - 1, i] * CDy + C[Ny - 3, i] * (CAy + CDy) + \
        C[Ny - 2, i + 1] * CDx + C[Ny - 2, i - 1] * (CAx + CDx) - C[Ny - 2, i] \
        * (CAx + CAy + 2 * CDx + 2 * CDy)
        
    # Inner portion of the domain
    for i in range(2, Ny - 2):
    
        for j in range(2, Nx - 2):
            
            sp[i ,j] = C[i + 1, j] * CDy + C[i - 1, j] * (CDy + 2 * CAy) - \
            C[i - 2, j] * (CAy / 2) + C[i, j + 1] * CDx + C[i, j - 1] * (CDx + \
            2 * CAx) - C[i, j - 2] * (CAx / 2) - C[i, j] * (2 * CDx + 2 * CDy \
            + 1.5 * CAy + 1.5 * CAx)
        
    return sp
        
# ==============================================================================
# Diffusion matrix assembly with Dirichlet Boundary conditions matching the 
# analytical solution at the boundary
# ==============================================================================
    
def Mat_ass(Sx, Sy, Nx, Ny):
    
    # Nodes with just one element in matrix
    N1 = (2 *Nx + 2 * Ny) - 4
    
    # Nodes with five elements in matrix 
    N2 = (Nx * Ny) - N1
    
    # Matrix elements corresponding to main diagonal Nx * Ny
    LHS_Dii = np.concatenate((np.ones(Nx + 1), np.ones(Nx - 1)), axis = 0)
    LHS_Iii = 
    LHS_Jii = 
    
    # Matrix elements correspnding to nodes (col + 1) and (col - 1)
    LHS_Di1 = 
    LHS_Ii1 = 
    LHS_Ji1 =
    LHS_Ji0 =
    
    # Matrix elements corresponding to (row + 1) and (row - 1)
    LHS_D1i = 
    LHS_I1i = 
    LHS_J1i =
    LHS_J0i =
    
    # Defining matrix as a sparse coordinate system matrix - empty lists
    LHS_d = np.zeros(N1 + 5* N2)
    LHS_i = np.zeros(N1 + 5* N2)
    LHS_j = np.zeros(N1 + 5* N2)
    
    # Filling main diagonal of data 
    LHS_d[0:N1] = np.ones()
    
    
    # Constructing coordinate matrix
    LHS = scsp.coo_matrix((Nx * Ny, Nx * Ny))    
    
    return LHS