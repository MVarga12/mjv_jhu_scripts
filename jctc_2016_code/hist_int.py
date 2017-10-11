#!/usr/bin/python

from __future__ import division
import os

print 'Before running this script, run integ.e/a.out on the lined up and corrected histogram to obtain integration of the full histogram'

#get integration from integ.e
i = input('Integration value: ')

hin = [line.strip() for line in open('histogramcorr.dat','r')]
hin = [float(hin) for hin in hin]
hinint = [x/i for x in hin]
histintfile = open('histogramcorrint.dat','w')
histintfile.writelines('%s \n' % t for t in hinint)

# get region for integration to obtain t he 'rate'

fin = hin[-60:]
finfile = open('productregion.dat','w')
finfile.writelines('%s \n' % t for t in fin)
