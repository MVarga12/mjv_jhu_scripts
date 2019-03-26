#!/usr/bin/env python2
##
# Created: 15-Jun-2018
# Modified: 16-Oct-2018
# Created by: Matthew Varga
#
# ## Purpose:
# Takes a trajectory, or single frame, its associated tpr file, and user provided
# atom pairs and returns a vector between the two atoms as well as:
#   1. The angle between the atom-pair vector and the z-axis, if only one pair is provided, or
#   2. The angle between the two atom-pair vectors, if two atom pairs are provided
# In the case of a trajectory, it then calculates the distribution of angles over the trajectory
# and computes the PDF. It fits the PDF to a function for nonharmonic angle constraints,
# as given by Gromacs (in section 4.1.3 in the 2016 release manual)
# The force constant, reference angle, and arbitrary constant for the constraint is returned.
# See (https://drive.google.com/file/d/1oDaAw1lzCjUYbac5tvmguXktIhurZxkG/view?usp=sharing) for
# more information
#
# ## Troubleshooting:
#   - If you get an error that reads, approximately, "Your tpx version is XXX, which this parser
#     does not support, yet.", you need to update MDAnalysis to the most recent version
#   - If being used on MARCC, or any other environment in which you cannot install python modules
#     globally and if you get errors which ultimately end up being due to MDAnalysis being out of date,
#     you will need to use virtualenv (pip install --user virtualenv) to set up a virtual environment
#     in which local packages are preferred to global packages (see script create_md_env.sh)
##

# without this, floating point division doesn't work properly
from __future__ import division

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import stats
from scipy.optimize import curve_fit
import argparse
import MDAnalysis  # for reading in and parsing trajectories
import numpy as np

# Use LaTeX for the labels, because the math looks a lot prettier
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


class Vector(object):
    # class to hold vectors created by atom pairs
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.mag = np.sqrt(x*x + y*y + z*z)

    # from an array of arrays like MDAnalysis.Universe gives you
    @classmethod
    def from_array(cls, arr):
        return cls(arr[0][0], arr[0][1], arr[0][2])

    # Normalizes the vector to unity
    def normalize(self):
        self.x = self.x / self.mag
        self.y = self.y / self.mag
        self.z = self.z / self.mag

    # Gets the cross product of the vector with another
    def cross(self, vec):
        u0 = self.x * vec.z - self.z * vec.y
        u1 = self.z * vec.x - self.x * vec.z
        u2 = self.x * vec.y - self.y * vec.x
        vecOut = Vector(u0, u1, u2)
        vecOut.normalize()
        return vecOut

    # returns the dot product of this Vector and another
    def dot(self, vec):
        return (self.x * vec.x) + (self.y * vec.y) + (self.z * vec.z)

    # returns the angle between this Vector and another
    def dot_theta(self, vec):
        theta = self.dot(vec)/(self.mag * vec.mag)

        # Precision tolerances
        if np.abs(theta - 1) < 1E-12:
            theta = 1
        if np.abs(-1 - theta) < 1E-12:
            theta = -1

        if self.mag == 0 or vec.mag == 0:
            return 0
        else:
            return np.arccos(theta)

    # Displays the vector in format [x, y, z]
    def display(self):
        print "[" + str(round(self.x, 5)) + ", " + \
            str(round(self.y, 5)) + ", " + str(round(self.z, 5)) + "]"


# Get the vector between the two atoms, while also accounting for possible periodic boundary conditions
def get_vector(pairFrame, pair):
    vector = Vector.from_array(pair[1] - pair[0])

    # correct for PBC
    corrX = vector.x - \
        (pairFrame.dimensions[0]*round(vector.x/pairFrame.dimensions[0]))
    corrY = vector.y - \
        (pairFrame.dimensions[1]*round(vector.y/pairFrame.dimensions[1]))
    corrZ = vector.z - \
        (pairFrame.dimensions[2]*round(vector.z/pairFrame.dimensions[2]))
    return Vector(corrX, corrY, corrZ)


# Gets the angle between two vectors, while also taking into account the sign
# of the angle.
def get_angle(vec1, vec2):
    ang = vec1.dot_theta(vec2)

    # Determine the sign by projecting the vectors onto the y-z plane
    # and determining the sign of the x component of the cross product
    tmpVec1 = Vector(0, vec1.y, vec1.z)
    tmpVec2 = Vector(0, vec2.y, vec2.z)
    if (tmpVec1.cross(tmpVec2).x > 0):
        return -ang
    else:
        return ang


# Non-harmonic angle constraint function, from Gromacs Manual, section 4.1.3
# E(d) = (1/2) * k (d - d_0)^2
def angle_fit(x, k, theta_0, c):
    return k*(1-np.cos(((np.pi/180)*(x-theta_0)))) + c


# Just displays the results, the two vectors and the angle between them
def display_results(pairList, arg, vec1, vec2=None):
    print "\nRESULTS:"
    print "Vector 1: atoms (" + pairList[0][0] + \
        ") and (" + pairList[0][1] + "):\n\t",
    vec1.display()

    if vec2 is not None:
        print "Vector 2: atoms (" + pairList[1][1] + \
            ") and (" + pairList[1][1] + "):\n\t",
        vec2.display()
    else:
        print "Vector 2: z-Axis vector"

    print "Angle (radians): ", ang
    print "Angle (degrees): ", ang * (180 / np.pi)


# Plots the histogram and probability distribution function of the angles collected along the trajectory.
# See documentation at beginning of script for specific information
def plot_histogram(angles, fname):
    # Get angles and calculate the number of bins
    minAng, maxAng = min(angles), max(angles)
    numBins = int(np.floor(2*(len(angles)**(1/3))))  # Rice rule
    print "Minimum angle:", round(minAng, 2)
    print "Maximum angle:", round(maxAng, 2)

    # Create and save the histogram of angles
    print "Creating histogram with", numBins
    hist = plt.hist(angles, density=True, bins=numBins, alpha=0.8)

    # Get a more manageable number of tick labels
    xt = plt.xticks()[0]
    xMin, xMax = min(xt), max(xt)
    lnspc = np.linspace(xMin, xMax, len(angles))

    # Save the histogram
    plt.title(r"{} angle histogram".format(fname.replace('_', ' ')))
    plt.xlabel(r"\text{Angle (degrees)}")
    plt.ylabel(r"\text{Counts}")

    # Create output file names
    outFName1, outFName2 = "", ""
    if not fname:
        outFName1 = r"angle_histogram.png"
        outFName2 = r"angle_invert.png"
    else:
        outFName1 = fname + r".png"
        outFName2 = fname + r"_invert.png"

    plt.savefig(outFName1, bbox_inches="tight", dpi=300)
    plt.close()

    # Get the "energy", E(d) \sim - kT*ln(PDF)
    # get bin centers instead of bin edges
    binCenters = 0.5*(hist[1][1:] + hist[1][:-1])
    invertHist = -np.log(hist[0])*0.008314462*310  # numpy's ln is numpy.log
    histX = hist[0]
    histY = hist[1]

    # Delete zeros, so that the "inverted" histogram is fit properly
    delete = []
    for i in range(0, len(invertHist)):
        if (np.abs(invertHist[i]) == np.inf):
            delete.append(i)
    delete.sort(reverse=True)  # delete indices from largest to smallest
    for i in delete:
        invertHist = np.delete(invertHist, i)
        binCenters = np.delete(binCenters, i)
        histX = np.delete(histX, i)
        histY = np.delete(histY, i)

    # Plot PDF data (histogram inverted)
    plt.scatter(binCenters, invertHist, marker='x', s=7, label=r"\text{data}")

    # Fit the PDF data to the function for angle constriants, given by Gromacs (section 4.1.3 in manual)
    k, theta_0, c = curve_fit(angle_fit, binCenters,
                              invertHist, bounds=(0, [2000, 360., np.inf]))[0]
    print k, theta_0, c, hist[1][np.argmax(hist[0])]
    plt.plot(lnspc, angle_fit(lnspc, k, theta_0, c), c='r', ls='--', lw=1,
             label=r"fit, k = {0:.3f}, $\theta_{{0}}$ = {1:.3f}, c = {2:.3f}".format(k, theta_0, c))

    # Save the PDF figure
    plt.legend()
    plt.xlabel("Angle ($^{\circ}$)")
    plt.ylabel("$-\log(PDF)*kT$")
    plt.savefig(outFName2, bbox_inches="tight", dpi=300)
    plt.close()

    # Make a polar histogram of the angles
    # Plots the histogram in polar coordinates. Makes it easier to visualize
    polAng = np.array(binCenters)*(np.pi/180)
    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location('N')
    bars = ax.bar(polAng, histX, width=(
        np.max(polAng) - np.min(polAng))/numBins, bottom=0)
    ax.set_xticklabels(["$0$", "$\\frac{\pi}{2}$", "$\\frac{\pi}{2}$", "$\\frac{3\pi}{4}$",
                        "$\pi$", "$-\\frac{3\pi}{4}$", "$-\\frac{\pi}{2}$", "$-\\frac{\pi}{4}$"], fontsize=18)
    ax.set_rticks([max(hist[0])/4, max(hist[0])/2, (3*max(hist[0]))/4])

    # Make the plot prettier by taking colours from a colourmap
    col = polAng - np.min(polAng)
    col != max(col)
    ax.set_rlabel_position(-135)
    ax.grid(linestyle='--')
    for r, bar in zip(col, bars):
        bar.set_facecolor(plt.cm.gnuplot2(r))
        bar.set_alpha(0.8)
    plt.savefig(outFName1[:-4] + "_polar.png", dpi=300, bbox_inches='tight')
    plt.close()

    # save the histogram and angle/constraint plot data
    with open(outFName1 + '.dat', 'w') as wout:
        for i in range(0, len(hist[0])):
            wout.write("%7.3f    %7.3f\n" % (hist[1][i], hist[0][i]))

    with open(outFName2 + '.dat', 'w') as wout:
        for i in range(0, len(binCenters)):
            wout.write("%7.3f     %7.3f\n" % (binCenters[i], invertHist[i]))

    # Print the results of the fit
    print "k*(1-cos(theta - theta_0)) + c)"
    print "k: %7.3f, theta_0: %7.3f, c: %7.3f" % (k, theta_0, c)


def read_frame(frame, system, sele1, sele2, frames, angles, wOut):
    # Read the frame and save angle to file
    vec = get_vector(system, [system.select_atoms(
        sele1).positions, system.select_atoms(sele2).positions])
    ang = get_angle(vec, Vector(0, 0, 1))
    angles.append(ang*(180/np.pi))
    frames.append(frame)
    wOut.write("%8i \t %8.4f \t %8.4f \t %8.4f \t %5.2f\n" %
               (frame, vec.x, vec.y, vec.z, ang*(180/np.pi)))


# set up argument parser and groups
parser = argparse.ArgumentParser(description='Calculate the vector and angle between one or two sets of atoms, using \
        an trajectory or GRO format frame and their associated TPR file.')
reqdArgs = parser.add_argument_group('required arguments')
frameArgs = parser.add_argument_group('frame arguments')
trajArgs = parser.add_argument_group('trajectory arguments')

# set up arguments
reqdArgs.add_argument('-p', '--tpr', action='store', dest='tpr',
                      help="TPR file for trajectory or single frame")
reqdArgs.add_argument('-pr', '--pairs', action='store', dest='pairs', type=int,
                      help="Number of pairs (1 or 2) to get get vectors and angles from")
parser.add_argument('-o', action='store', dest='fname',
                    type=str, help='Prefix of output filenames.')
trajArgs.add_argument('-n', '--frame', action='store', dest='frame', type=int,
                      help="Frame number (timestep) from which you want to extract the angle(s) and vector(s)")
trajArgs.add_argument('-t', '--traj', action='store', dest='traj',
                      help="XTC or TRR format trajectory file")
trajArgs.add_argument('-in', action='store', dest='intr', type=int,
                      help='Data point interval if calculating vectors and angles along the entire trajectory.')
frameArgs.add_argument('-g', '--gro', action='store',
                       dest='gro',  help="GRO single frame file")
parser.add_argument('-pf', '--pair_file', action='store', dest='pairFile',
                    help='File containing atom pairs, in VMD syntax, with one atom per line.')

args = parser.parse_args()

# quit conditions, for if information is missing at runtime
if args.tpr is None or args.tpr[-3:] != 'tpr':
    print "Error, please include a valid TPR file"
    exit()

if args.gro is not None and args.traj is not None and args.frame is not None:
    print "Error, please only use an XTC file with frame number or GRO file alone, not both."
    exit()

if args.traj is not None:
    if (args.intr is None and args.frame is None):
        print "Error, trajectory detected but no frame or data point interval is provided."
        exit(1)
    if (args.intr is not None and args.frame is not None):
        print "Error, only provide single frame number or data point interval."
        exit(1)

# Get the atom pair(s) from the file, or from user input at runtime
# MDAnalysis uses similar syntax to VMD, with a few exceptions (bynum instead of index)
print "Using VMD syntax\n"
rawPairs = []
if args.pairFile is not None:
    print "Reading from pairs file..."
    tmpPairs = [line.strip() for line in open(args.pairFile, 'r')]
    rawPairs.append([tmpPairs[0], tmpPairs[1]])
    if len(tmpPairs) > 2:
        rawPairs.append([tmpPairs[2], tmpPairs[3]])
elif args.pairs == 1:
    rawPairs.append([raw_input("Enter the first atom: "),
                     raw_input("Enter the second atom: ")])
elif args.pairs == 2:
    rawPairs.append([raw_input("Enter the first atom of the first pair: "), raw_input(
        "Enter the second atom of the first pair: ")])
    rawPairs.append([raw_input("Enter the first atom of the second pair: "), raw_input(
        "Enter the second atom of the second pair: ")])
else:
    print "Please provide one or two valid atom pair(s) in VMD syntax."
    exit()

# Create systems and get desired frame
if args.gro is not None:
    print "Reading frame from", args.gro
    pairFrame = MDAnalysis.Universe(args.tpr, args.gro)
else:
    if args.intr is None:
        # IF LOOKING FOR ONE FRAME
        print "Reading trajectory from", args.traj,
        pairFrame = MDAnalysis.Universe(args.tpr, args.traj)
        pairFrame.trajectory[args.frame]
        print "with " + str(len(pairFrame.trajectory)) + " frames."
    else:
        # IF LOOKING ALONG ENTIRE TRAJECTORY
        system = MDAnalysis.Universe(args.tpr, args.traj)

        # Set up output files, with names dependent on user provided information
        # e.g. filename and number of pairs
        # If filename not provided, name based on the name of the trajectory provided
        if (not args.fname) and (args.pairs == 1):
            outFName = "vecs_angs_" + \
                args.traj.split("/")[-1].split(".")[0] + ".dat"
        elif (not args.fname) and (args.pairs == 2):
            outFName = "vecs_angs_" + \
                args.traj.split("/")[-1].split(".")[0] + "_pair1.dat"
        elif args.fname and args.pairs == 1:
            outFName = args.fname + ".dat"
        elif args.fname and args.pairs == 2:
            outFName = args.fname + "_pair1.dat"
        else:
            print "ERROR IN CREATING FILENAMES"
            exit(1)

        wOut = open(outFName, 'w')
        angles, frames = [], []

        # if only one pair is provided, find the angle between the vector between them and the z-axis
        if args.pairs == 1:
            wOut.write("(%s) \t (%s)\n" % (rawPairs[0][0], rawPairs[0][1]))
            wOut.write("%8s \t %8s \t %8s \t %8s \t %5s\n" %
                       ("ts", "x", "y", "z", "ang"))  # print output header

            for frame in range(0, len(system.trajectory), args.intr):
                # Skip the first 500ns
                if frame < 5000:
                    continue

                system.trajectory[frame]  # ffwd to frame number
                read_frame(
                    frame, system, rawPairs[0][0], rawPairs[0][1], frames, angles, wOut)

            # Plot the histogram/PDF
            plot_histogram(angles, args.fname)

            # Plot angle of the vectors over the course of the trajectory
            plt.plot(frames, angles, lw=1, c='black')
            plt.xlabel("Frame")
            plt.ylabel("Angle (degrees)")
            plt.savefig(args.fname + "_angle_timeseries.png",
                        bbox_inches="tight", dpi=300)

            exit()  # we're done, exit so we don't try to go further
        else:
            # If two pairs of atoms are provided, get the angles between the vectors created by them.
            angles2 = []
            if not args.fname:
                outFName = "vecs_angs_" + \
                    args.traj.split("/")[-1].split(".")[0] + "_pair2.dat"
            else:
                outFName = args.fname + "_pair2.dat"

            wOut2 = open(outFName, "w")
            wOut2.write("(%s) \t (%s)\n" % (rawPairs[1][0], rawPairs[1][1]))
            wOut2.write("%8s \t %8s \t %8s \t %8s \t %5s\n" %
                        ("ts", "x", "y", "z", "ang"))  # print header

            for frame in range(0, len(system.trajectory), args.intr):
                # Skip the first 500ns
                if frame < 5000:
                    continue

                system.trajectory[frame]  # ffwd to frame number
                read_frame(
                    frame, system, rawPairs[0][0], rawPairs[0][1], frames, angles, wOut)

                # PAIR 2
                dummy = []  # Only want to append the frame number to the frame list once
                read_frame(
                    frame, system, rawPairs[1][0], rawPairs[1][1], dummy, angles2, wOut2)

            allAngles = angles + angles2
            plot_histogram(allAngles, args.fname)

            # Plot angle over time
            plt.plot(frames, angles, lw=1, c='black', label="Pair 1")
            plt.plot(frames, angles2, lw=1, c='red', label="Pair 2")
            plt.legend()
            plt.title("Angle with z-axis over trajectory")
            plt.xlabel("Frame")
            plt.ylabel("Angle (degrees)")
            plt.savefig(args.fname + "_angle_timeseries.png",
                        bbox_inches="tight", dpi=300)

            exit()  # we're done, exit so we don't try to go further

# IF ONLY USING ONE FRAME
# Get pair coordinates from the desired frame
if args.pairs == 1:
    # if only one pair is provided, find the angle between the vector between them and the z-axis
    vec1 = get_vector(pairFrame, [pairFrame.select_atoms(
        rawPairs[0][0]).positions, pairFrame.select_atoms(rawPairs[0][1]).positions])
    ang = get_angle(vec1, Vector(0, 0, 1))
    display_results(rawPairs, ang, vec1, None)
else:
    # If two pairs are provided, get the angle between their vectors
    vec1 = get_vector(pairFrame, [pairFrame.select_atoms(
        rawPairs[0][0]).positions, pairFrame.select_atoms(rawPairs[0][1]).positions])
    vec2 = get_vector(pairFrame, [pairFrame.select_atoms(
        rawPairs[1][0]).positions, pairFrame.select_atoms(rawPairs[1][1]).positions])
    ang = get_angle(vec1, vec2)
    display_results(rawPairs, ang, vec1, vec2)
