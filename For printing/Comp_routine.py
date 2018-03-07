#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Programa para comparar los resultados de esquemas
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana
Febrero de 2018
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import T2D 

# =============================================================================
# Defining parameters for plotting and running
# =============================================================================

t0 = 1.
tf = 5. 

# =============================================================================
# Running with 21 nodes (dx = dy = 0.25 m)
# =============================================================================

dt = 0.1
time = np.linspace(t0, tf, 40)
ee01 = T2D.AB2D(21, 21, dt)
eft01_1 = T2D.FT2D(21, 21, dt, 1.)
ef01_05 = T2D.FT2D(21, 21, dt, 0.5)

style.use('ggplot')
plt.figure(1)
#plt.subplot(1, 1, 1)
plt.semilogy(time, ee01, '-r', label = 'Explicit')
plt.semilogy(time, eft01_1, '-b', label = 'Fractional timestep implicit')
plt.semilogy(time, ef01_05, '-g', label = 'Fractional timestep CN')
plt.title('Norm 2 Error evolution')
plt.legend()
plt.ylim([1e-5, 1e7])
plt.savefig('Error01.pdf')

dt = 0.01
time = np.linspace(t0, tf, 400)
ee001 = T2D.AB2D(21, 21, dt)
eft001_1 = T2D.FT2D(21, 21, dt, 1.)
ef001_05 = T2D.FT2D(21, 21, dt, 0.5)

plt.figure(2)
#plt.subplot(1, 1, 1)
plt.semilogy(time, ee001, '-r', label = 'Explicit')
plt.semilogy(time, eft001_1, '-b', label = 'Fractional timestep implicit')
plt.semilogy(time, ef001_05, '-g', label = 'Fractional timestep CN')
plt.title('Norm 2 Error evolution')
plt.legend()
plt.ylim([1e-5, 1e-1])
plt.savefig('Error001.pdf')

dt = 0.001
time = np.linspace(t0, tf, 4000)
ee0001 = T2D.AB2D(21, 21, dt)
eft0001_1 = T2D.FT2D(21, 21, dt, 1.)
ef0001_05 = T2D.FT2D(21, 21, dt, 0.5)

plt.figure(3)
#plt.subplot(1, 1, 1)
plt.semilogy(time, ee0001, '-r', label = 'Explicit')
plt.semilogy(time, eft0001_1, '-b', label = 'Fractional timestep implicit')
plt.semilogy(time, ef0001_05, '-g', label = 'Fractional timestep CN')
plt.title('Norm 2 Error evolution')
plt.legend()
plt.ylim([1e-5, 1e-1])
plt.savefig('Error0001.pdf')

dt = 0.0001
time = np.linspace(t0, tf, 40000)
ee00001 = T2D.AB2D(21, 21, dt)
eft00001_1 = T2D.FT2D(21, 21, dt, 1.)
ef00001_05 = T2D.FT2D(21, 21, dt, 0.5)

plt.figure(4)
#plt.subplot(1, 1, 1)
plt.semilogy(time, ee00001, '-r', label = 'Explicit')
plt.semilogy(time, eft00001_1, '-b', label = 'Fractional timestep implicit')
plt.semilogy(time, ef00001_05, '-g', label = 'Fractional timestep CN')
plt.title('Norm 2 Error evolution')
plt.legend()
plt.ylim([1e-5, 1e-1])
plt.savefig('Error00001.pdf')



disc_tiempo = np.array([0.1, 0.01, 0.001, 0.0001])
e_exp = np.array([np.linalg.norm(ee01), np.linalg.norm(ee001), \
                  np.linalg.norm(ee0001), np.linalg.norm(ee00001)])
e_imp = np.array([np.linalg.norm(eft01_1), np.linalg.norm(eft001_1), \
                  np.linalg.norm(eft0001_1), np.linalg.norm(eft00001_1)])
e_cn = np.array([np.linalg.norm(ef01_05), np.linalg.norm(ef001_05), \
                  np.linalg.norm(ef0001_05), np.linalg.norm(ef00001_05)])

plt.figure(5)
plt.loglog(disc_tiempo, e_exp, '-r', label='Explicit')
plt.loglog(disc_tiempo, e_imp, '-b', label='Fractional timestep implicit')
plt.loglog(disc_tiempo, e_cn, '-g', label='Fractions timestep CN')
plt.legend()
plt.savefig('Convergencia.pdf')

