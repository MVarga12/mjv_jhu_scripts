#!/usr/bin/env python2
# Purpose: plots the two proteins' center of mass and the z-axis-constrained atom to box edge z-axis

import MDAnalysis
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
import sys

## FIGURE TEXT STUFF ##
# matplotlib.rcParams['lines.linewidth'] = 1
# matplotlib.rcParams['font.family'] = 'Garuda'
# matplotlib.rcParams['text.usetex'] = True

## GET DATA ##
system = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
pro1 = system.residues[0:214]
pro2 = system.residues[214:428]
pro1Atom = pro1.atoms[336]
pro2Atom = pro2.atoms[336]

frames, pro1ComDist, pro2ComDist, pro1AtomDist, pro2AtomDist = [], [], [], [], []
for ts in system.trajectory:
    frames.append(ts.frame)
    pro1ComDist.append(0.1*np.abs(pro1.center_of_mass()[2] - system.dimensions[2]))
    pro2ComDist.append(0.1*np.abs(pro2.center_of_mass()[2] - system.dimensions[2]))
    pro1AtomDist.append(0.1*np.abs(pro1Atom.position[2] - system.dimensions[2]))
    pro2AtomDist.append(0.1*np.abs(pro2Atom.position[2] - system.dimensions[2]))

## PLOT ##
fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10,4))
ax[0][0].plot(frames,pro1ComDist, label="Pro1 COM\n(%7.3f $\pm$ %5.3f nm)" % (np.mean(pro1ComDist), np.std(pro1ComDist)), c='pink')
ax[1][0].plot(frames,pro2ComDist, label="Pro2 COM\n(%7.3f $\pm$ %5.3f nm)" % (np.mean(pro2ComDist), np.std(pro2ComDist)), c='lightblue')
ax[0][1].plot(frames,pro1AtomDist, label="Pro1 155LEU\n(%7.3f $\pm$ %5.3f nm)" % (np.mean(pro1AtomDist), np.std(pro1AtomDist)), c='#800517')
ax[1][1].plot(frames,pro2AtomDist, label="Pro2 155LEU\n(%7.3f $\pm$ %5.3f nm)" % (np.mean(pro2AtomDist), np.std(pro2AtomDist)), c='darkblue')

## FORMATTING ##
fig.tight_layout()
plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=4, fontsize='small', fancybox=True, shadow=True)
fig.text(0.5,0,"Frame",ha='center')
fig.text(0,0.6,"Distance to box edge\n(z-axis, nm)",ha='center',rotation='vertical')
# sns.despine(top=True, right=True)
fig.savefig("test.png", dpi=300, bbox_inches='tight')
