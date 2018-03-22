#!/usr/local/bin/python
'''
MJV: March 2018
quick script to read in data files of coordinates and plot them
    in 3D space as scatterplot.
cycles between colours, so that each protein has a unique colour
'''
import matplotlib.pyplot as plt
import argparse
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(description='plot protein crds from dat files')
parser.add_argument('input files', metavar='I', type=str, nargs='+', help='input file for plotting')
args = parser.parse_args()

filenames = []
for k in args.__dict__:
    if args.__dict__[k] is not None:
        for j in args.__dict__[k]:
            filenames.append(str(j))

color=iter(plt.cm.rainbow(np.linspace(0,1,len(filenames))))

#import data
cfig = plt.figure()
clat = cfig.add_subplot(111, projection='3d')

for i in range(0, len(filenames)):
    filein = [line.split() for line in open(filenames[i], 'r')]

    x, y, z = [], [], []
    for i in range(1, len(filein)):
            for j in range(0, len(filein[i])):
                filein[i][j] = float(filein[i][j])
            x.append(filein[i][0])
            y.append(filein[i][1])
            z.append(filein[i][2])

    clat.scatter(x, y, z, c=next(color), label = i)

clat.set_xlabel("x")
clat.set_ylabel("y")
clat.set_zlabel("z")

plt.show()
