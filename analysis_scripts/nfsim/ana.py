#!/usr/bin/env python2
# script to take columns from a formatted data file 
# and return the averages of several runs
#import sys

#filename = str(sys.argv[1])

d = [line.split() for line in open("things_bound.dat","r")]

sum1=0
sum2=0
sum3=0
for i in range(0,3):
    sum1+=float(d[i+1][2])
    sum2+=float(d[i+1][3])
    sum3+=float(d[i+1][4])

sum1 = sum1/3
sum2 = sum2/3
sum3 = sum3/3

print sum1, sum2, sum3
