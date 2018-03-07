##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Ejemplo de 4 figuras en una para implementar en el caso de transporte bidimensio
#nal en la ecuacion de transporte. 
#
#@author: toni
#"""
#
## Four axes, returned as a 2-d array
#f, axarr = plt.subplots(2, 2)
#axarr[0, 0].plot(x, y)
#axarr[0, 0].set_title('Axis [0,0]')
#axarr[0, 1].scatter(x, y)
#axarr[0, 1].set_title('Axis [0,1]')
#axarr[1, 0].plot(x, y ** 2)
#axarr[1, 0].set_title('Axis [1,0]')
#axarr[1, 1].scatter(x, y ** 2)
#axarr[1, 1].set_title('Axis [1,1]')
## Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
#
#
## Just a figure and one subplot
#f, ax = plt.subplots()
#ax.plot(x, y)
#ax.set_title('Simple plot')
#
## Two subplots, the axes array is 1-d
#f, axarr = plt.subplots(2, sharex=True)
#axarr[0].plot(x, y)
#axarr[0].set_title('Sharing X axis')
#axarr[1].scatter(x, y)
#
## Two subplots, unpack the axes array immediately
#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#ax1.plot(x, y)
#ax1.set_title('Sharing Y axis')
#ax2.scatter(x, y)
#
## Three subplots sharing both x/y axes
#f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#ax1.plot(x, y)
#ax1.set_title('Sharing both axes')
#ax2.scatter(x, y)
#ax3.scatter(x, 2 * y ** 2 - 1, color='r')
## Fine-tune figure; make subplots close to each other and hide x ticks for
## all but bottom plot.
#f.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#
## row and column sharing
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#ax1.plot(x, y)
#ax1.set_title('Sharing x per column, y per row')
#ax2.scatter(x, y)
#ax3.scatter(x, 2 * y ** 2 - 1, color='r')
#ax4.plot(x, 2 * y ** 2 - 1, color='r')
#
#
#"Ejemplo colorbar0"
#https://matplotlib.org/examples/pylab_examples/colorbar_tick_labelling_demo.html
#
#import matplotlib.pyplot as plt
#import numpy as np
#from matplotlib import cm
#from numpy.random import randn
#
## Make plot with vertical (default) colorbar
#fig, ax = plt.subplots()
#
#data = np.clip(randn(250, 250), -1, 1)
#
#cax = ax.imshow(data, interpolation='nearest', cmap=cm.coolwarm)
#ax.set_title('Gaussian noise with vertical colorbar')
#
## Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar
#
## Make plot with horizontal colorbar
#fig, ax = plt.subplots()
#
#data = np.clip(randn(250, 250), -1, 1)
#
#cax = ax.imshow(data, interpolation='nearest', cmap=cm.afmhot)
#ax.set_title('Gaussian noise with horizontal colorbar')
#
#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
#cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
#
#plt.show()
#
#
## Librerias para graficar que fueron importadas en algún momento 
#import matplotlib.pyplot as plt
#from matplotlib import style
#from matplotlib import cm
#import matplotlib.animation as animation


import matplotlib.pyplot as plt
import numpy as np

plt.ion(); plt.figure(1);
for k in range(10):
    plt.clf(); plt.subplot(121);
    plt.contourf(np.random.randn(10,10)); plt.colorbar();
    plt.subplot(122,polar=True)
    plt.contourf(np.random.randn(10,10)); plt.colorbar();
    plt.draw();
    plt.pause(0.2)