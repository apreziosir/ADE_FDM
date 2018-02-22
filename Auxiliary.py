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