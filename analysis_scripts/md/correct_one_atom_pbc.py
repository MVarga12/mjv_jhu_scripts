#!/usr/bin/env python2
import argparse
import MDAnalysis
from MDAnalysis.analysis import distances

parser = argparse.ArgumentParser(
    description="Script to extract one frame from a trajectory.")
reqdArgs = parser.add_argument_group("required arguments")
reqdArgs.add_argument('-p', '--tpr', action='store', dest='tpr', type=str,
                      help="TPR file corresponding to the trajectory (str)", required=True)
reqdArgs.add_argument('-g', '--gro', action='store', dest='gro', type=str,
                      help="GRO file (str)", required=True)
args = parser.parse_args()

system = MDAnalysis.Universe(args.tpr, args.gro)
atom = system.select_atoms("bynum " + raw_input("Which atom (index): "))

corr_pos = []
for i in range(0, 3):
    corr_pos.append(atom.positions[0][i] - (atom.dimensions[i]
                                            * round(atom.positions[0][i] / atom.dimensions[i])))

print corr_pos
