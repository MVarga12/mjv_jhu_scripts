#!/usr/bin/env python2
##
# Created: 13-Sep-2018 10:06:27 AM EDT
# Modified: 13-Sep-2018 01:54:02 PM EDT
# Created by: Matthew Varga
# Purpose: Plot 2D free energy surface from plumed sum_hills output
##

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import numpy as np
import sys

# Get data
print "Reading data..."
data = [line.split() for line in open(sys.argv[1], "r")]
data2 = [x for x in data if not x == []]  # strip out headers

d1, d2, free, dd1, dd2 = [], [], [], [], []
for elem in data2[9:]:
    d1.append(float(elem[0]))
    d2.append(float(elem[1]))
    free.append(float(elem[2]))
    dd1.append(float(elem[3]))
    dd2.append(float(elem[4]))

X = np.linspace(min(d1), max(d1), 1318)
Y = np.linspace(min(d2), max(d2), 1322)

print "Creating data grid. This may take a while..."
D1, D2 = np.meshgrid(X, Y)
ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)

levels = np.arange(np.min(free), np.max(free), 10)

contour = plt.contour(D1, D2, ENER, colors='k', linewidths=0.5, levels=levels)
contourf = plt.contourf(D1, D2, ENER, cmap=cm.Spectral, levels=levels)
plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel("Distance D1 (nm)")
plt.ylabel("Distance D2 (nm)")
plt.savefig("fes2d.png", dpi=300, bbox_inches='tight')
plt.close()

# 3D plots
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(D1, D2, ENER, cmap=cm.Spectral)
ax.set_zlim(np.min(ENER) - (np.max(ENER) - np.min(ENER))/10)
ax.contour(D1, D2, ENER, 10, linewidths=1, levels = levels, linestyles='solid', offset=ax.get_zlim()[0], cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='y', offset=np.min(D2), cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='x', offset=np.min(D1), cmap=cm.Spectral)
ax.set_zlabel("Free Energy (kJ/mol)")
ax.set_ylabel("Distance D2 (nm)")
ax.set_xlabel("Distance D1 (nm)")
ax.view_init(20, 45)
plt.draw()
plt.savefig("3dfes1.png", dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(D1, D2, ENER, cmap=cm.Spectral)
ax.set_zlim(np.min(ENER) - (np.max(ENER) - np.min(ENER))/10)
ax.contour(D1, D2, ENER, 10, linewidths=1, levels = levels, linestyles='solid', offset=ax.get_zlim()[0], cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='y', offset=np.max(D2), cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='x', offset=np.min(D1), cmap=cm.Spectral)
ax.set_zlabel("Free Energy (kJ/mol)")
ax.set_ylabel("Distance D2 (nm)")
ax.set_xlabel("Distance D1 (nm)")
ax.view_init(20, -45)
plt.draw()
plt.savefig("3dfes2.png", dpi=300, bbox_inches='tight')
