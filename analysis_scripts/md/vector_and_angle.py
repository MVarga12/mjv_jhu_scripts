#!/usr/bin/env python2
##
# Created: 15-Jun-2018
# Modified: 20-Jun-2018
# Created by: Matthew Varga
# Purpose: Translation of vector_and_angle.m into python, using mdanalysis
# Takes a trajectory, or single frame, its associated tpr file, and user provided
# atom pairs and returns a vector between the two atoms as well as:
#   1) the angle between the atom-pair vector and the z-axis, if only one pair is provided, or
#   2) the angle between the two atom-pair vectors, if two atom pairs are provided
#
# Troubleshooting:
#   - If you get an error that reads, approximately, "Your tpx version is XXX, which this parser
#     does not support, yet.", you need to update MDAnalysis to the most recent version
##

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats
import argparse
import MDAnalysis
import numpy as np


class Vector(object):
    # class to hold vectors
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.mag = np.sqrt(x*x + y*y + z*z)

    # from an array of arrays like Universe gives you
    @classmethod
    def from_array(cls, arr):
        return cls(arr[0][0], arr[0][1], arr[0][2])

    # returns the dot product of this Vector and another
    def dot(self, vec):
        return (self.x * vec.x) + (self.y * vec.y) + (self.z * vec.z)

    # returns the angle between this Vector and another
    def dot_theta(self, vec):
        theta = self.dot(vec)/(self.mag * vec.mag)

        if np.abs(theta - 1) < 1E-12:
            theta = 1
        if np.abs(-1 - theta) < 1E-12:
            theta = -1

        if self.mag == 0 or vec.mag == 0:
            return 0
        else:
            return np.arccos(theta)

    def display(self):
        print "[" + str(round(self.x, 5)) + ", " + \
            str(round(self.y, 5)) + ", " + str(round(self.z, 5)) + "]"


def get_vector(pairFrame, pair):
    # Get the vector between the two atoms
    vector = Vector.from_array(pair[1] - pair[0])

    # correct for PBC
    corrX = vector.x - \
        (pairFrame.dimensions[0]*round(vector.x/pairFrame.dimensions[0]))
    corrY = vector.y - \
        (pairFrame.dimensions[1]*round(vector.y/pairFrame.dimensions[1]))
    corrZ = vector.z - \
        (pairFrame.dimensions[2]*round(vector.z/pairFrame.dimensions[2]))
    return Vector(corrX, corrY, corrZ)


def get_angle(vector1, vector2):
    return vector1.dot_theta(vector2)


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


def plot_histogram(angles):
    # Bin Size
    binSize = 0.2

    # Get angles and calculate the number of bins
    minAng, maxAng = min(angles), max(angles)
    #numBins = int(np.floor((maxAng - minAng)/binSize))
    numBins = int(np.floor(2*(len(angles)**(1/3))))  # Rice rule
    print len(angles), numBins
    print "Minimum angle:", round(minAng, 2)
    print "Maximum angle:", round(maxAng, 2)

    # Create and save the histogram
    print "Creating histogram with", numBins, "bins, using bin size", binSize
    hist = plt.hist(angles, normed=True, bins=numBins)

    xt = plt.xticks()[0]
    xMin, xMax = min(xt), max(xt)
    lnspc = np.linspace(xMin, xMax, len(angles))

    # Gaussian
    mean, stdev = stats.norm.fit(angles)
    gaussPDF = stats.norm.pdf(lnspc, mean, stdev)
    plt.plot(lnspc, gaussPDF, label="Norm")

    # Gamma
    ag, bg, cg = stats.gamma.fit(angles)
    gammaPDF = stats.gamma.pdf(lnspc, ag, bg, cg)
    plt.plot(lnspc, gammaPDF, label="Gamma")

    # Beta
    ab, bb, cb, db = stats. beta.fit(angles)
    betaPDF = stats.beta.pdf(lnspc, ab, bb, cb, db)
    plt.plot(lnspc, betaPDF, label="Beta")

    plt.legend()
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Counts")
    plt.savefig("angle_histogram.png", bbox_inches="tight", dpi=300)

    # with open("angle_histogram.dat", 'w') as wOut:
    #    for point in hist:
    #        wOut.write("%7.3f \t %7.3f\n" % (point[0], point[1]))

    #print "Most likely angle:", hist[1][np.argmax(hist[0])]
    print "Most likely angle (Gaussian):", lnspc[np.argmax(gaussPDF)]
    print "Most likely angle (Gamma):", lnspc[np.argmax(gammaPDF)]
    print "Most likely angle (Beta):", lnspc[np.argmax(betaPDF)]


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

# quit conditions
# TODO: add checks for xtc and gro file types
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

# Get the atom pair(s)
print "Using VMD syntax\n"
rawPairs = []
if args.pairFile is not None:
    print "Reading from pairs file..."
    tmpPairs = [line.strip() for line in open(args.pairFile, 'r')]
    rawPairs.append([tmpPairs[0], tmpPairs[1]])
    if len(rawPairs) == 2:
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

        angles = []
        with open("vecs_angs_" + args.traj.split("/")[-1].split(".")[0] + ".dat", 'w') as wOut:
            print "Writing vectors and angles to vecs_angs_" + \
                args.traj.split("/")[-1].split(".")[0] + ".dat."
            if args.pairs == 1:
                wOut.write("(%s) \t (%s)\n" % (rawPairs[0][0], rawPairs[0][1]))
                wOut.write("%8s \t %8s \t %8s \t %5s\n" %
                           ("x", "y", "z", "ang"))  # print header
                # if only one pair is provided, find the angle between the vector between them and the z-axis
                for frame in range(0, len(system.trajectory), args.intr):
                    system.trajectory[frame]  # ffwd to frame number
                    vec = get_vector(system, [system.select_atoms(
                        rawPairs[0][0]).positions, system.select_atoms(rawPairs[0][1]).positions])
                    ang = get_angle(vec, Vector(0, 0, 1))
                    angles.append(ang*(180/np.pi))
                    wOut.write("%8.4f \t %8.4f \t %8.4f \t %5.2f\n" %
                               (vec.x, vec.y, vec.z, ang*(180/np.pi)))
                plot_histogram(angles)
                exit()  # we're done, exit so we don't try to go further
            else:
                # print header
                wOut.write("(%s) \t (%s) \t\t (%s \t %s" % (
                    rawPairs[0][0], rawPairs[0][1], rawPairs[1][0], rawPairs[1][1]))
                # else, find the angle between the vectors created by each pair
                for frame in range(0, len(system.trajectory), args.intr):
                    system.trajectory[frame]  # ffwd to frame number
                    vec1 = get_vector(system, [system.select_atoms(
                        rawPairs[0][0]).positions, system.select_atoms(rawPairs[0][1]).positions])
                    vec2 = get_vector(system, [system.select_atoms(
                        rawPairs[1][0]).positions, system.select_atoms(rawPairs[1][1]).positions])
                    ang = get_angle(vec1, vec2)
                    angles.append(ang*(180/np.pi))
                    wOut.write(
                        "%8.4f \t %8.4f \t %8.4f \t %8.4f \t % 8.4f \t %8.4f \t %5.2f\n" % (vec1.x, vec1.y, vec1.z, vec2.x, vec2.x, vec2.x, ang * (180/np.pi)))
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
    vec1 = get_vector(pairFrame, [pairFrame.select_atoms(
        rawPairs[0][0]).positions, pairFrame.select_atoms(rawPairs[0][1]).positions])
    vec2 = get_vector(pairFrame, [pairFrame.select_atoms(
        rawPairs[1][0]).positions, pairFrame.select_atoms(rawPairs[1][1]).positions])
    ang = get_angle(vec1, vec2)
    display_results(rawPairs, ang, vec1, vec2)
