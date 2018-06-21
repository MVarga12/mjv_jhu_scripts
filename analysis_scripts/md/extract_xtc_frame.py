#!/usr/bin/env python2
# #
# Created:  7-May-2018 01:38:58 PM EDT
# Modified:  7-May-2018 02:51:17 PM EDT
# Created by: Matthew Varga
# Purpose: quick script to use MDAnalysis to extract one frame from an XTC trajectory
# runtime arguments are tpr_file xtc_file target_frame
##

import MDAnalysis
import os.path
import sys


def print_frame(system, inFile, frame, segmentChoice):
    system.trajectory[frame]
    segmentOut = system.select_atoms(segmentChoice)
    wout = MDAnalysis.Writer("frame_" + str(frame) + "_" + inFile, segmentOut.n_atoms)
    wout.write(segmentOut)


# get topology file
inTop = sys.argv[1]

# check if the file exists and has the right extension
if (inTop[-3:] != "tpr"):
    print "Error, not the correct extension."
    exit()

if not os.access(inTop, os.R_OK):
    print "Error, " + inTop + " doesn't exist."
    exit()

# get the input trajectory name
inFile = sys.argv[2]

# check if the file exists and has the right extension
if (inFile[-3:] != "xtc"):
    print "Error, not the correct extension."
    exit()

if not os.access(inFile, os.R_OK):
    print "Error, ", inFile, " doesn't exist."
    exit()

# get the target frame
frame = int(sys.argv[3])
if not isinstance(frame, (int, long)):
    print "Error, the provided frame is not an integer."
    exit()

# does the user want the trajectory to contain all atoms or just a portion
# would like to add in options here to output what segments are available
part = int(raw_input("Get all atoms [1] or specific segments [2]: "))

system = MDAnalysis.Universe(inTop, inFile)
if part == 1:
    print_frame(system, inFile, frame, "all")
elif part == 2:
    print_frame(system, inFile, frame, raw_input("Portion to select (VMD syntax): "))
else:
    print "Error, not a correct segment choice."
    exit()
