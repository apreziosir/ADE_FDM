#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auxiliary routines for the FDM FTCS 2D transport equation
@author: Antonio Preziosi-Ribero
CFD Pontificia Universidad Javeriana
"""

import time
import numpy as np

# ==============================================================================
# Printing stability crteria calculated. Remember what happens with Sx and Sy 
# in the explicit case, and take into account that none of the CFL can be higher
# than 1.0
# ==============================================================================

def show_sc(Sx, Sy, CFLx, CFLy):
    
    print('Stability criteria for the 2D transport problem given: \n')
    print('================================')
    print('Sx\t', 'Sy\t', 'CFLx\t', 'CFLy')
    print('================================')
    print('{:.5}'.format(Sx), '\t', '{:.5}'.format(Sy), '\t', 
          '{:.5}'.format(np.abs(CFLx)), '\t', '{:.5}'.format(np.abs(CFLy)))
    print('================================')
    
    time.sleep(2)

    return 