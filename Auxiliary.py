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

