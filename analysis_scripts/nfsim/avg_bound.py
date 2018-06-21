#!/usr/bin/env python2
import sys

filename = str(sys.argv[1])

d = [line.split() for line in open(filename,"r")]
del d[0]
for i in range(0,len(d)):
	for j in range(0,len(d[i])):
		d[i][j] = float(d[i][j])

for i in range(0,len(d)):
	cltc = d[i][3]/d[i][1]
	epn1 = d[i][4]/d[i][2]
        l = 1-(d[i][5]/1000)
	print "%.3f %.3f %.3f" % (cltc, epn1,l)
