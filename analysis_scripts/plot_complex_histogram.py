#!/usr/bin/env python2
# Created 06/06/19, MJV
# Plots the complex component histogram of the final frame
# ONLY WORKS FOR HOMOGENOUS COMPLEXES FROM HOMOGENOUS SIMULATIONS

import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

from matplotlib import rc
font = {'family' : 'sans-serif',
        'size' : 20 }
rc('font', **font)

parser = argparse.ArgumentParser(description="Plots the histogram of complex sizes.")
reqdArgs = parser.add_argument_group('required arguments')
optArgs = parser.add_argument_group('optional arguments')

# set up arguments
reqdArgs.add_argument('-f', action='append', dest='file', type=str, help='Complex histogram file, output from NERDSS as a file named ComplexHistogram_Np$NUMPROS.dat')
optArgs.add_argument('-x', '--xlabel', dest='xlabel', type=str, help='X-axis label', default='Species per complex')
optArgs.add_argument('-y', '--ylabel', dest='ylabel', type=str, help='Y-axis label', default='Counts')
optArgs.add_argument('-t', '--title', dest='title', type=str, help='Title of the plot', default=None)
args = parser.parse_args()

if len(args.file) == 1:
    rawData = [line.split() for line in open(args.file[0])]
    data = []
    header = "iter:"
    for elem in reversed(rawData):
        if (elem[0] != header):
            data.append([float(elem[0]), elem[1], float(elem[2])])
        else:
            break

    data.sort(key=lambda x: x[2]) # sort by the number of species in each complex
    data = zip(*data) # transpose the array
    plt.bar(data[2], data[0], width=0.2, alpha=0.8)
    plt.xticks(np.arange(1, data[2] + 1)) # make sure the x axis ticks are only integers
else:
    rawData = []
    for fileName in args.file:
        rawData.append([line.split() for line in open(fileName)])
    data = []
    maxIndex = 0
    for frame in rawData:
        tmpData = []
        for elem in reversed(frame):
            if elem[0] != "iter:":
                tmpData.append([float(elem[2]), float(elem[0])])
            else:
                break
        tmpData.sort(key=lambda x: x[0])
        tmpData = map(list, zip(*tmpData))
        if max(tmpData[0]) > maxIndex:
            maxIndex = max(tmpData[0])
        data.append(pd.DataFrame(tmpData[1], index=tmpData[0]))

    for frame in data:
        frame = frame.reindex(np.arange(1, maxIndex + 1))
        frame = frame.fillna(0)

    # Create a data frame containing the average values and SEM of each
    # complex size
    dfCat = pd.concat(data).groupby(level=0)
    avg = dfCat.mean()
    sem = dfCat.sem()
    comb = avg.join(sem, lsuffix="-mn", rsuffix="-sem").fillna(0)

    # create the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(comb.index.tolist(), comb['0-mn'], width=0.3, alpha=0.8, yerr=comb['0-sem'], capsize=4)
    plt.xticks(np.arange(1, comb.index.tolist()[-1] + 1)) # make sure the x axis ticks are only integers


# plot the data
plt.xlabel(args.xlabel)
plt.ylabel(args.ylabel)
if args.title:
    plt.title(args.title)
plt.savefig('complex_histogram.png', dpi=300, bbox_inches='tight')
