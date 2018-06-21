#!/usr/bin/env python2
# #
# Created:  7-May-2018
# Modified:  20-June-2018
# Created by: Matthew Varga
# Purpose: quick script to use MDAnalysis to extract one frame from an XTC trajectory
# runtime arguments are tpr_file xtc_file target_frame
##

import MDAnalysis
import argparse
import os.path
import sys


def print_frame(system, traj, frame, segChoice):
    # iterate to correct frame (yes, it isn't set equal to anything on purpose)
    system.trajectory[frame]

    # select desired atoms to write
    try:
        segOut = system.select_atoms(segChoice)
    except:
        print "Error, problem with atom selection syntax."
        exit(1)

    # prepare writer, using file name with the form "frame_<frame number>_<traj name>.gro"
    wout = MDAnalysis.Writer("frame_" + str(frame) +
                             "_" + traj.split("/")[-1].split(".")[0] + ".gro", reindex=False)
    wout.write(segOut)


parser = argparse.ArgumentParser(
    description="Script to extract one frame from a trajectory.")
reqdArgs = parser.add_argument_group("required arguments")
reqdArgs.add_argument('-p', '--tpr', action='store', dest='tpr', type=str,
                      help="TPR file corresponding to the trajectory (str)", required=True)
reqdArgs.add_argument('-t', '--traj', action='store', dest='traj', type=str,
                      help="Trajectory containing the frame for extraction (str)", required=True)
reqdArgs.add_argument('-n', '--frame', action='store', dest='frame',
                      type=long, help="Number of the frame for extraction. Remember that VMD starts lists from index 1 *shudder*", required=True)
args = parser.parse_args()

# check if the file exists and has the right extension
if (args.tpr[-3:] != "tpr"):
    print "Error, not the correct extension."
    exit()

if not os.access(args.tpr, os.R_OK):
    print "Error, " + args.tpr + " doesn't exist."
    exit()

# check if the file exists and has the right extension
if (args.traj[-3:] != "xtc" and args.traj[-3:] != "trr"):
    print "Error, not the correct extension."
    exit()

if not os.access(args.traj, os.R_OK):
    print "Error, ", args.traj, " doesn't exist."
    exit()

if not isinstance(args.frame, (int, long)):
    print "Error, please provide a frame number."
    exit()

# does the user want the trajectory to contain all atoms or just a portion
# would like to add in options here to output what segments are available
system = MDAnalysis.Universe(args.tpr, args.traj)
print_frame(system, args.traj, args.frame, raw_input(
    "Portion to select: "))
