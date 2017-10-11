#!/usr/bin/python

from __future__ import division
from subprocess import call
import collections
import os

lines = [line.split() for line in open('histogramcorr.dat','r')]
op = [float(line[0]) for line in lines]
#print op
corrtot = [float(line[1]) for line in lines]
#print corrtot
intin = raw_input('Enter full histogram integral:' )
intin = float(intin)

# INTEGRATION

def trapezoid(f,a,b,n):
      h = (b-a)/n
      s = 0.0
      s += f[a]/2.0
      for i in range(1,n):
            s += f[a+1]
      s += f[b]/2.0
      return s * h

#integral = trapezoid(corrtot,0,49,49)
#print trapezoid(lambda x:x**2,5,10,100)

hcorr = [s/intin for s in corrtot]
histcorr = open('histogramcorr_int.dat','w')
histcorr.writelines('%s %s \n' % t for t in zip(op,hcorr))
histcorr.close()

print len(corrtot)
print corrtot[38]
# with full bins (11) hcorr should be hcorr[76:]

prod = hcorr[37:]
opx = op[37:]
prodhist = open('productregion.dat','w')
prodhist.writelines('%s %s \n' % t for t in zip(opx,prod))
prodhist.close()
print op[37]
