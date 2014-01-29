#Use numpy and scipy
import math
import numpy as np 
import scipy as sy 
from matplotlib import pyplot

#Read data from first file into two arrays
#could also use array = np.loadtext()
x1=[]
y1=[]
for line in open('/home/becca/astr596/data/sn2011by-hst+lick.flm'):
	columns = line.split()
	if len(columns) == 2:
		if float(columns[0]) <= 5715: # truncation
			x1.append(float(columns[0]))
			y1.append(float(columns[1]))

#Read data from second file into two arrays

x2 = []
y2 = []
for line in open('/home/becca/astr596/data/sn2011fe-visit3-hst.flm'):
	columns = line.split()
	if len(columns) == 2:
		x2.append(float(columns[0]))
		y2.append(float(columns[1]))

#printed data check that it read in right
#print x1
#print y1
#print x2
#print y2

#Plot arrays one by the other
pyplot.plot(x1, y1, color='c', linewidth=2)
pyplot.plot(x2, y2, color='m', linewidth=2)
#pyplot.axis([float('1e-17'), float('1e-12'), 1200, 6000]) #define range for x and y
pyplot.yscale('log')
pyplot.show()