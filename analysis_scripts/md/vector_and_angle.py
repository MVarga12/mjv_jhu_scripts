#!/usr/local/bin/python
##
# Created: 15-Jun-2018 12:08:57 PM EDT
# Modified: 18-Jun-2018 04:23:58 PM EDT
# Created by: Matthew Varga
# Purpose: Translation of vector_and_angle.m into python, using mdanalysis
# Takes a trajectory, or single frame, its associated tpr file, and user provided
# atom pairs and returns a vector between the two atoms as well as:
#   1) the angle between the atom-pair vector and the z-axis, if only one pair is provided, or
#   2) the angle between the two atom-pair vectors, if two atom pairs are provided
##

import MDAnalysis
import argparse
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
    def dot(cls, vec):
        return (cls.x * vec.x) + (cls.y * vec.y) + (cls.z * vec.z)

    # returns the angle between this Vector and another
    def dot_theta(cls, vec):
        theta = cls.dot(vec)/(cls.mag * vec.mag)

        if np.abs(theta - 1) < 1E-12:
            theta = 1
        if np.abs(-1 - theta) < 1E-12:
            theta = -1

        if cls.mag == 0 or vec.mag == 0:
            return 0
        else:
            return np.arccos(theta)

    def display(cls):
        print "[" + str(cls.x) + ", " + str(cls.y) + ", " + str(cls.z), "]"


def get_vector(pair):
    # Get the vector between the two atoms
    vector = Vector.from_array(pair[1] - pair[0])

    # correct for PBC
    corrX = vector.x - (system.dimensions[0]*round(vector.x/system.dimensions[0]))
    corrY = vector.y - (system.dimensions[1]*round(vector.y/system.dimensions[1]))
    corrZ = vector.z - (system.dimensions[2]*round(vector.z/system.dimensions[2]))
    return Vector(corrX, corrY, corrZ)


def get_angle(system, vector1, vector2):
    print ""


# set up argument parser and groups
parser = argparse.ArgumentParser(description='Calculate the vector and angle between one or two sets of atoms, using \
        an XTC format trajectory or GRO format frame and their associated TPR file.')
reqdArgs = parser.add_argument_group('required arguments')
frameArgs = parser.add_argument_group('frame arguments')
trajArgs = parser.add_argument_group('trajectory arguments')

# set up arguments
reqdArgs.add_argument('--tpr', action='store', dest='tpr', help="TPR file for trajectory or single frame")
reqdArgs.add_argument('--pairs', action='store', dest='pairs', type=int, help="Number of pairs (1 or 2) to get get vectors and angles from")
trajArgs.add_argument('--frame', action='store', dest='frame', type=int, default=-1, help="Frame number (timestep) from which you want to extract the angle(s) and vector(s)")
trajArgs.add_argument('--traj', action='store', dest='traj',  help="XTC format trajectory file")
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
    print "Reading from frame..."
    system = MDAnalysis.Universe(args.tpr, args.gro)
    pairFrame = system
else:
    print "Reading from trajectory..."
    system = MDAnalysis.Universe(args.tpr, args.traj)
    pairFrame = system.trajectory[args.frame]

# Get pair coordinates from the desired frame
pairCrds = []
for pair in rawPairs:
    pairCrds.append([pairFrame.select_atoms(pair[0]).positions, pairFrame.select_atoms(pair[1]).positions])
