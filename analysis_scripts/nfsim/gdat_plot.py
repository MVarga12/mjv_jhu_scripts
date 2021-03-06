#!/usr/bin/env python2
import sys
import subprocess

filename = sys.argv[1]  # gdat file to be analyzed
target = sys.argv[2]    # target specie

d = [line.split() for line in open(filename, "r")]
del d[0][0]  # delete that stupid pound sign

ans = raw_input("Alter column (y/n)?")
if ans == 'y':
    cons = raw_input("Divide column by what symmetry? (integer) ")
    cons = int(cons)
else:
    cons = 1

val = 0
for i in range(0, len(d[0])):
    if d[0][i] == target:
        val = i

time, targlist, out = [], [], []
for i in range(1, len(d)):
    time.append(d[i][0])

for i in range(1, len(d)):
    targlist.append(d[i][val])

for i in range(0, len(time)):
    out.append([float(time[i]), float(targlist[i])/cons])

f = open('tmp', 'w')
for i in range(0, len(out)):
    f.write("%f \t %f \n" % (out[i][0], out[i][1]))

subprocess.call("xmgrace tmp &", shell=True)
