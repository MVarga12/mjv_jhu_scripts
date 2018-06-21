#!/usr/bin/env python2
##
# Created: 15-Jun-2018 12:08:57 PM EDT
# Modified: 18-Jun-2018 04:23:58 PM EDT
# Created by: Matthew Varga
# Purpose: Translation of vector_and_angle.m into python, using mdanalysis
# Takes a trajectory, or single frame, its associated tpr file, and user provided
# atom pairs and returns a vector between the two atoms as well as:
#   1) the angle between the atom-pair vector and the z-axis, if only one pair is provided, or
#   2) the angle between the two atom-pair vectors, if two atom pairs are provided
#
# TODO: Change rawPairs/whatever to a class (stops it from being hardcoded as a list)
##

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
        print "[" + str(round(self.x, 5)) + ", " + str(round(self.y, 5)) + ", " + str(round(self.z, 5)) + "]"


def get_vector(pairFrame, pair):
    # Get the vector between the two atoms
    vector = Vector.from_array(pair[1] - pair[0])

    # correct for PBC
    corrX = vector.x - (pairFrame.dimensions[0]*round(vector.x/pairFrame.dimensions[0]))
    corrY = vector.y - (pairFrame.dimensions[1]*round(vector.y/pairFrame.dimensions[1]))
    corrZ = vector.z - (pairFrame.dimensions[2]*round(vector.z/pairFrame.dimensions[2]))
    return Vector(corrX, corrY, corrZ)


def get_angle(vector1, vector2):
    return vector1.dot_theta(vector2)


def display_results(pairList, arg, vec1, vec2 = None):
    print "\nFound angle between the two vectors."
    print "Vector 1: atoms (" + pairList[0][0] + ") and (" + pairList[0][1] + "):\n\t",
    vec1.display()

    if vec2 is not None:
        print "Vector 2: atoms (" + pairList[1][1] + ") and (" + pairList[1][1] + "):\n\t",
        vec2.display()
    else:
        print "Vector 2: z-Axis vector"

    print "Angle (radians): ", ang
    print "Angle (degrees): ", ang * (180/np.pi)

# set up argument parser and groups
parser = argparse.ArgumentParser(description='Calculate the vector and angle between one or two sets of atoms, using \
        an trajectory or GRO format frame and their associated TPR file.')
reqdArgs = parser.add_argument_group('required arguments')
frameArgs = parser.add_argument_group('frame arguments')
trajArgs = parser.add_argument_group('trajectory arguments')

# set up arguments
reqdArgs.add_argument('--tpr', action='store', dest='tpr', help="TPR file for trajectory or single frame")
reqdArgs.add_argument('--pairs', action='store', dest='pairs', type=int, help="Number of pairs (1 or 2) to get get vectors and angles from")
trajArgs.add_argument('--frame', action='store', dest='frame', type=int, default=-1, help="Frame number (timestep) from which you want to extract the angle(s) and vector(s)")
trajArgs.add_argument('--traj', action='store', dest='traj',  help="XTC or TRR format trajectory file")
frameArgs.add_argument('--gro', action='store', dest='gro',  help="GRO single frame file")
parser.add_argument('--pair_file', action='store', dest='pairFile', help='File containing atom pairs, in VMD syntax, with one atom per line.')

args = parser.parse_args()

# quit conditions
# TODO: add checks for xtc and gro file types
if args.tpr is None or args.tpr[-3:] != 'tpr':
    print "Error, please include a valid TPR file"
    exit()

if (args.gro is None) and (args.traj is None or args.frame == -1):
    print "Error, please include a valid GRO file or a valid XTC file with a frame number."
    exit()

if args.gro is not None and args.traj is not None and args.frame != -1:
    print "Error, please only use an XTC file with frame number or GRO file alone, not both."
    exit()

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
    rawPairs.append([raw_input("Enter the first atom: "), raw_input("Enter the second atom: ")])
elif args.pairs == 2:
    rawPairs.append([raw_input("Enter the first atom of the first pair: "), raw_input("Enter the second atom of the first pair: ")])
    rawPairs.append([raw_input("Enter the first atom of the second pair: "), raw_input("Enter the second atom of the second pair: ")])
else:
    print "Please provide one or two valid atom pair(s) in VMD syntax."
    exit()

# Create systems and get desired frame
if args.gro is not None:
    print "Reading frame from", args.gro
    pairFrame = MDAnalysis.Universe(args.tpr, args.gro)
else:
    print "Reading trajectory from", args.traj,
    pairFrame = MDAnalysis.Universe(args.tpr, args.traj)
    print "with " + str(len(pairFrame.trajectory)) + " frames."

# Get pair coordinates from the desired frame
if args.pairs == 1:
    # if only one pair is provided, find the angle between the vector between them and the z-axis
    vec1 = get_vector(pairFrame, [pairFrame.select_atoms(rawPairs[0][0]).positions, pairFrame.select_atoms(rawPairs[0][1]).positions])
    ang = get_angle(vec1, Vector(0, 0, 1))
    display_results(rawPairs, ang, vec1, None)
else:
    vec1 = get_vector(pairFrame, [pairFrame.select_atoms(rawPairs[0][0]).positions, pairFrame.select_atoms(rawPairs[0][1]).positions])
    vec2 = get_vector(pairFrame, [pairFrame.select_atoms(rawPairs[1][0]).positions, pairFrame.select_atoms(rawPairs[1][1]).positions])
    an = get_angle(vec1, vec2)
    display_results(rawPairs, ang, vec1, vec2)
