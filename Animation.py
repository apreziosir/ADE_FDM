#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:37:25 2018

@author: apreziosir
"""

"""
=========================
Simple animation examples
=========================

This example contains two animations. The first is a random walk plot. The
second is an image animation.
"""

#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#
#
#def update_line(num, data, line):
#    line.set_data(data[..., :num])
#    return line,
#
#fig1 = plt.figure()
#
#data = np.random.rand(2, 25)
#l, = plt.plot([], [], 'r-')
#plt.xlim(0, 1)
#plt.ylim(0, 1)
#plt.xlabel('x')
#plt.title('test')
#line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#                                   interval=50, blit=True)
#
## To save the animation, use the command: line_ani.save('lines.mp4')
#
#fig2 = plt.figure()
#
#x = np.arange(-9, 10)
#y = np.arange(-9, 10).reshape(-1, 1)
#base = np.hypot(x, y)
#ims = []
#for add in np.arange(15):
#    ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))
#
#im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
#                                   blit=True)
## To save this second animation with some metadata, use the following command:
## im_ani.save('im.mp4', metadata={'artist':'Guido'})
#
#plt.show()
#
#

#
#"""
#=================
#An animated image
#=================
#
#This example demonstrates how to animate an image.
#"""
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#
#fig = plt.figure()
#
#
#def f(x, y):
#    return np.sin(x) + np.cos(y)
#
#x = np.linspace(0, 2 * np.pi, 120)
#y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
#
#im = plt.imshow(f(x, y), animated=True)
#
#
#def updatefig(*args):
#    global x, y
#    x += np.pi / 15.
#    y += np.pi / 20.
#    im.set_array(f(x, y))
#    return im,
#
#ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
#plt.show()





"""
============
3D animation
============

#A simple example of an animated plot... In 3D!
#"""
#import numpy as np
#import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d.axes3d as p3
#import matplotlib.animation as animation
#
#
#def Gen_RandLine(length, dims=2):
#    """
#    Create a line using a random walk algorithm
#
#    length is the number of points for the line.
#    dims is the number of dimensions the line has.
#    """
#    lineData = np.empty((dims, length))
#    lineData[:, 0] = np.random.rand(dims)
#    for index in range(1, length):
#        # scaling the random numbers by 0.1 so
#        # movement is small compared to position.
#        # subtraction by 0.5 is to change the range to [-0.5, 0.5]
#        # to allow a line to move backwards.
#        step = ((np.random.rand(dims) - 0.5) * 0.1)
#        lineData[:, index] = lineData[:, index - 1] + step
#
#    return lineData
#
#
#def update_lines(num, dataLines, lines):
#    for line, data in zip(lines, dataLines):
#        # NOTE: there is no .set_data() for 3 dim data...
#        line.set_data(data[0:2, :num])
#        line.set_3d_properties(data[2, :num])
#    return lines
#
## Attaching 3D axis to the figure
#fig = plt.figure()
#ax = p3.Axes3D(fig)
#
## Fifty lines of random 3-D lines
#data = [Gen_RandLine(25, 3) for index in range(50)]
#
## Creating fifty line objects.
## NOTE: Can't pass empty arrays into 3d version of plot()
#lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
#
## Setting the axes properties
#ax.set_xlim3d([0.0, 1.0])
#ax.set_xlabel('X')
#
#ax.set_ylim3d([0.0, 1.0])
#ax.set_ylabel('Y')
#
#ax.set_zlim3d([0.0, 1.0])
#ax.set_zlabel('Z')
#
#ax.set_title('3D Test')
#
## Creating the Animation object
#line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),
#                                   interval=50, blit=False)
#
#plt.show()




"""
=====
Decay
=====

This example showcases a sinusoidal decay animation.
"""
#
#
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#
#
#def data_gen(t=0):
#    cnt = 0
#    while cnt < 1000:
#        cnt += 1
#        t += 0.1
#        yield t, np.sin(2*np.pi*t) * np.exp(-t/10.)
#
#
#def init():
#    ax.set_ylim(-1.1, 1.1)
#    ax.set_xlim(0, 10)
#    del xdata[:]
#    del ydata[:]
#    line.set_data(xdata, ydata)
#    return line,
#
#fig, ax = plt.subplots()
#line, = ax.plot([], [], lw=2)
#ax.grid()
#xdata, ydata = [], []
#
#
#def run(data):
#    # update the data
#    t, y = data
#    xdata.append(t)
#    ydata.append(y)
#    xmin, xmax = ax.get_xlim()
#
#    if t >= xmax:
#        ax.set_xlim(xmin, 2*xmax)
#        ax.figure.canvas.draw()
#    line.set_data(xdata, ydata)
#
#    return line,
#
#ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,
#                              repeat=False, init_func=init)
#plt.show()





import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 6*np.pi, 100)
y = np.sin(x)

# You probably won't need this if you're embedding things in a tkinter plot...
plt.ion()

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'r-') # Returns a tuple of line objects, thus the comma

for phase in np.linspace(0, 10*np.pi, 500):
    line1.set_ydata(np.sin(x + 3))
    fig.canvas.draw()