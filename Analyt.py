#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytical solution to 2D transport equation
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana 
"""

import numpy as np 

# ==============================================================================
# Analytical solution to de two dimensional transport equation in a rectangular
# domain with constant velocity fields
# ==============================================================================

def difuana(M, A, Dx, Dy, u, v, x0, y0, x, y, t):
    
    temp1 = A * 4 * np.pi * (t) * np.sqrt(Dx * Dy)
    
    temp2 = ((x - x0 - u * t) ** 2) / (4 * Dx * t)
    
    temp3 = ((y - y0 - v * t) ** 2) / (4 * Dy * t)
    
    C = (M / temp1) * np.exp(-temp2 - temp3)
    
    return C 