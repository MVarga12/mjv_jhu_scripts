#!/usr/bin/env python2
# this reads a .species file from nfsim and looks at each species
# for the CLTC, EPN1, and L
# it then finds the number of CLTC and EPN1 in each species
# and reports the size of each aggregate, iff it is located in the same
# species as L, i.e. if the aggregate is on the membrane

from __future__ import division
import sys
from itertools import groupby

filename = str(sys.argv[1])
out = str(sys.argv[2])

run = [line.split() for line in open(filename, "r")]
del run[0]
del run[0]

species = []
for i in range(0,len(run)):
    x = run[i][0].split(".")
    x.append(run[i][1])
    species.append(x)

vals = []
cltc = "CLTC"
lip = "L"

#extract number of cltcs in species, and if the species is on the membrane
for i in range(0,len(species)):
    cltc_num = 0
    lip_num = 0
    spe_num = 0
    for j in range(0,len(species[i])-1):
        if species[i][j][0:4] == cltc:
            cltc_num += 1
        if species[i][j][0] == lip:
            lip_num += 1
    if lip_num >= 1:
        lip_num = 1
    spe_num = species[i][len(species[i])-1]
    x = [spe_num,cltc_num,lip_num]
    vals.append(x)


# remove species not on the membrane
for i in range(0,len(vals)):
    if (vals[i][1] == 0 or vals[i][2] == 0):
        vals[i] = []

vals = [elem for elem in vals if elem != []]

cltc_hist = []
for i in range(0,len(vals)):
    for j in range(int(vals[i][0])):
        cltc_hist.append(vals[i][1])

cltc_hist.sort()
fin = [(key, len(list(group))) for key, group in groupby(cltc_hist)]

#tuple to list of lists of ints
fin2 = []
for i in range(0,len(fin)):
    x = int(fin[i][0])
    y = int(fin[i][1])
    tmp = [x,y]
    fin2.append(tmp)
    tmp = []

print species[len(species)-1]

#normalization
denom = 0
for i in range(0,len(fin2)):
    denom += fin2[i][1]

for i in range(0,len(fin2)):
    tmp = float(fin2[i][1])/float(denom) #normalize
    fin2[i].append(tmp)
    tmp = fin2[i][0] * fin2[i][1]       #total clat in agg size
    fin2[i].append(tmp)

#d = open(out,"w")
#for i in range(0,len(fin2)):
#    print >> d, '%i \t %i \t %f \t %i' % (fin2[i][0],fin2[i][1],fin2[i][2],fin2[i][3])
