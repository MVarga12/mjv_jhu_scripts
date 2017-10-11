
from __future__ import division
from subprocess import call
import collections
import os
from scipy import stats

# Definitions
endwin = []
win_num = 8
dict = {}
dict2 = {}
dictop = {}
dictfin = {}

for i in range(1,win_num):
      if i not in dict:
            dict[i] = list()
      if i not in dict2:
            dict2[i] = list()
      if i not in dictop:
            dictop[i] = list()

# Define order parameter dict for output of final histogram

for i in range(1, win_num):
      if i == 1:
            dictop[1] = range(200,360,10)
#      if i > 1:
#            dictop[i] = range(dictop[i-1][len(dictop[i-1])-1] - 5, dictop[i-1][len(dictop[i-1])-3] + 65, 5)
      if i == 2:
            dictop[i] = range(350,410,10)
      if i == 7:
            dictop[i] = range(600,900,10)
      if i == 3:
            dictop[i] = range(400,460,10)
      if i == 4:
            dictop[i] = range(450,510,10)
      if i == 5:
            dictop[i] = range(500,560,10)
      if i == 6:
            dictop[i] = range(550,610,10)
      print dictop[i]

# Create order parameter files from CHARMM outputted distances files. 
# NOTE: Order parameter in this instance is the fractional hydride-donor distance
# NOTE: If error that list 'aint' is too small, delete the most recently
# created order parameter file and run again
'''
for j in range(1,win_num):
      if not os.path.exists('op_'+str(j)):
            os.makedirs('op_'+str(j))

for j in range(1,win_num):
      for i in range(2,10000,1):
            try:
                  if not os.path.exists('op_'+str(j)+'/op_'+str(j)+'_'+str(i)+'.dat'):
                        x = [line.strip() for line in open('../window_'+str(j)+'/distances/disha_'+str(j)+'_'+str(i)+'.dat','r')]
                        y = [line.strip() for line in open('../window_'+str(j)+'/distances/dishd_'+str(j)+'_'+str(i)+'.dat','r')]
                        xfloat = [float(x) for x in x]
                        yfloat = [float(y) for y in y]
                        op = [y/(x+y) for x,y in zip(xfloat,yfloat)]
                        outfile = open('op_'+str(j)+'/op_'+str(j)+'_'+str(i)+'.dat','w')
                        outfile.writelines('%s \n' % u for u in op)
                        outfile.close()
            except IOError:
                  break
#                  m = i
#                  for i in range(m,10000,1):
#                        try:
#                              if not os.path.exists('op_'+str(j)+'/op_'+str(j)+'_'+str(i)+'.dat'):
#                                    x = [line.strip() for line in open('../window_'+str(j)+'/distances/disha_'+str(j)+'_'+str(i)+'.dat','r')]
#                                    y = [line.strip() for line in open('../window_'+str(j)+'/distances/dishd_'+str(j)+'_'+str(i)+'.dat','r')]
#                                    xfloat = [float(x) for x in x]
#                                    yfloat = [float(y) for y in y]
#                                    op = [y/(x+y) for x,y in zip(xfloat,yfloat)]
#                                    outfile = open('op_'+str(j)+'/op_'+str(j)+'_'+str(i)+'.dat','w')
#                                    outfile.writelines('%s \n' % u for u in op)
#                                    outfile.close()
#                        except IOError:
#                              m = m+1
#                              if m-i == 21:
#                                    break
# Read the endpoints of the previously created order parameter files into lists within
# dict2 for later processing
'''

for j in range(1,win_num):
      for i in range(2,10000,1):
            try:
                  a = [line.strip() for line in open('../new_bins/window_'+ str(j)+'/op_'+str(j)+'/op_'+str(j)+'_'+str(i)+'.dat','r')] 
                  aint = [float(b) for b in a]
#                 print len(aint)
                  dict[j].append(aint[250])
            except IOError:
                  break
for i in range (1,win_num):
      print i, len(dict[i])
# Make probability distributions of the histograms
# NOTE: Each window will have smaller windows (bins) that can be changed on the fly
# by changing the step of the range of i

def histogram(start,stop,win):
      for i in range(start,stop,10):
            m=0
            for j in range(1,len(dict[win])):
                  #print i/1000, dict[win][j], (i/1000) + 0.005
                  if (i/1000) <= dict[win][j] <= (i/1000)+0.010:
                        m = m + 1
            dict2[win].append(m)
      dict2[win] = [float(s/sum(dict2[win])) for s in dict2[win]]
      print 'Window '+str(win)+' completed!'

def histogram2(start,stop,win):
      for i in range(start,stop,25):
            m=0
            for j in range(1,len(dict[win])):
                  #print i/1000, dict[win][j], (i/1000) + 0.005
                  if (i/1000) <= dict[win][j] <= (i/1000)+0.025:
                        m = m + 1
            dict2[win].append(m)
      dict2[win] = [float(s/sum(dict2[win])) for s in dict2[win]]
      print 'Window '+str(win)+' completed!'

#print dict[1]
histogram(200,360,1)
histogram(350,410,2)
histogram(400,460,3)
histogram(450,510,4)
histogram(500,560,5)
histogram(550,610,6)
histogram(600,900,7)
#histogram(550,605,8)
#histogram(600,900,9)i

#print dict2[4]
#histogram(start,stop,win)

for i in range(1,win_num):
      hist = open('histogram_'+str(i)+'.dat','w')
      hist.writelines('%s %s \n' % t for t in zip(dictop[i],dict2[i]))
      hist.close()

# Create a 'corrected' histogram by matching up the endpoints of the histograms
# in the overlap regions by taking the ratio and multiplying the second
# histogram by that ratio

print dict2[2]
dict2[2][0] = (dict2[2][0] + dict2[2][1]) / 2
for i in range(2,5):
      dict2[2][i-1] = dict2[2][i]
print dict2[2]      
del dict2[2][3]
dict2[2] = [float(s/sum(dict2[2])) for s in dict2[2]]
print dict2[2]
print sum(dict2[2])
print dictop[2]
dictop[2][0] = 360
for i in range(2,5):
      dictop[2][i-1] = dictop[2][i]
del dictop[2][3]
print dictop[2]

print dict2[2]
dict2[2][0] = (dict2[2][0] + dict2[2][1]) / 2
for i in range(2,5):
      dict2[2][i-1] = dict2[2][i]
print dict2[2]      
del dict2[2][3]
dict2[2] = [float(s/sum(dict2[2])) for s in dict2[2]]
print dict2[2]
print sum(dict2[2])
print dictop[2]
dictop[2][0] = 370
for i in range(2,5):
      dictop[2][i-1] = dictop[2][i]
del dictop[2][3]
print dictop[2]

print "2:", dict2[2]
'''
print dict2[6]
dict2[6][1] = (dict2[6][1] + dict2[6][2]) / 2
for i in range(3,5):
      dict2[6][i-1] = dict2[6][i]
del dict2[6][4]
dict2[6] = [float(s/sum(dict2[6])) for s in dict2[6]]
print dict2[6]
print sum(dict2[6])
print dictop[6]
dictop[6][1] = 560
for i in range(3,5):
      dictop[6][i-1] = dictop[6][i]
del dictop[6][4]
print dictop[6]

for i in range(1,3):
      if i == 1:
            y = dict2[i]
            x = [x/1000 for x in dictop[i]]
            x[len(x)-1] = x[len(x)-1] + 0.005
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            #print dictop[i], x, x[0], x[len(x)-1]
            matchval2 = slope*x[len(x)-1] + intercept
            dict2[i][len(dict2[i])-1] = matchval2
#      elif i == 6:
#            continue             
      else:
            y = dict2[i]
            x = [x/1000 for x in dictop[i]]
            x[0] = x[0] + 0.005
            x[len(x)-1] = x[len(x)-1] + 0.005
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            matchval1 = slope*x[0] + intercept
            matchval2 = slope*x[len(x)-1] + intercept
            dict2[i][0] = matchval1
            dict2[i][len(dict2[i])-1] = matchval2 
      print i, slope, dict2[i][0], dict2[i][int(len(dict2[i])-1)]

'''
corrtot = dict2[1]
#for i in range(1,win_num-1):
#      print i
#      if i ==1:
#            corr = dict2[i][len(dict2[i])-1] / dict2[i+1][0]
#            new = [(s*corr) for s in dict2[i+1]]
#      elif i > 1:
#            print new[len(new)-1], dict2[i+1][0]
#            if new[len(new)-1]==0:
#                  corr = 2.0
#            if new[len(new)-1]==0 and dict2[i+1][0]==0:
#                  corr = 0.5 
#                  print "puppies"
#            else:
#                  corr = new[len(new)-1] / dict2[i+1][0]
#                  print corr
#            new = [(s*corr) for s in dict2[i+1]]
#      corrtot.extend(new)

for i in range(1,win_num-1):
      print i
      if i == 1:
            corr = dict2[i][len(dict2[i])-1] / dict2[i+1][0]
            print dict2[1][len(dict2[1])-1] , dict2[i+1][0]      
            new = [(s*corr) for s in dict2[i+1]]
            print i,corr
            del new[0]
#      elif i == 5:
#            new = dict2[i+1]
#            corr = new[len(new)-1] / dict2[i+1][0]
#            new = [(s*corr) for s in dict2[i+1]]
#            print corr, new      
      else:
            corr = new[len(new)-1] / dict2[i+1][0]
            print new[len(new)-1] , dict2[i+1][0]      
            new = [(s*corr) for s in dict2[i+1]]
            print i,corr
            del new[0]
      hist22 = open('histogram_'+str(i+1)+'_corr.dat','w')
      hist22.writelines('%s \n' % t for t in new)
      hist22.close()
      corrtot.extend(new)

op3 = dictop[1]
for i in range(2,win_num):
      del dictop[i][0]
      op3.extend(dictop[i])
      print op3
op3 = [x/1000 for x in op3]

final = open('histogramcorr.dat','w')
final.writelines('%s %s \n' % s for s in zip(op3,corrtot))
final.close()

print len(corrtot)
print len(dict2[7])
# INTEGRATION

#def trapezoid(f,a,b,n):
#      h = (b-a)/n
#      s = 0.0
#      s += f[a]/2.0
#      for i in range(1,n):
#            s += f[a+1]
#      s += f[b]/2.0
#      return s * h

#integral = trapezoid(corrtot,0,49,49)
#print trapezoid(lambda x:x**2,5,10,100)

#hcorr = [s/integral for s in corrtot]
#histcorr = open('histogramcorr_int.dat','w')
#histcorr.writelines('%s \n' % t for t in hcorr)
#histcorr.close()
'''
print trapezoid(hcorr,39,67,25)
'''
