#Use numpy and scipy
import math
import numpy as np 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#Read data from first file into two arrays
#could also use array = np.loadtext()
x1=[]
y1=[]
for line in open('/Users/Rebecca/astr596/data/sn2011by-hst+lick.flm'):
	columns = line.split()
	if len(columns) == 2:
		if float(columns[0]) <= 5715: # truncation
			x1.append(float(columns[0]))
			y1.append(float(columns[1]))

#Read data from second file into two arrays
#could also use array = np.loadtext()
x2 = []
y2 = []
for line in open('/Users/Rebecca/astr596/data/cfa/sn1995Y/sn1995Y-19950916.49-fast.flm'):
	columns = line.split()
	if len(columns) == 2:
		x2.append(float(columns[0]))
		y2.append(float(columns[1]))

def avg(*vals):
	return float(sum(vals))/float(len(vals))
#print avg(1,2,3,9.5) #test line for avg function

x3 = []
y3 = []
for x1_val in range(0,len(x1)):
	for x2_val in range (0,len(x2)):
		if x1[x1_val] == x2[x2_val]:
			x3.append(x1[x1_val])
			y3.append(avg(y1[x1_val],y2[x2_val]))
			x3_string = str(x3)
			y3_string = str(y3)
		x2_val += 1
	x1_val += 1

f3 = interp1d(x3, y3, kind='linear')
x4 = x3
y4 = f3(x3)

#--------------- the point at which my code becomes useful again
#Read averaged data to file -- change this to .flm or whatever
np.savetxt('averaged_spectra.txt', np.transpose([x3,y3]), fmt="%d %26.18e")

#Plot arrays one by the other
plt.subplot(111)
plt.plot(x1, y1, color='c', linewidth=1)
plt.plot(x2, y2, color='m', linewidth=1)
plt.plot(x3, y3, color='#ff6600', linewidth=1)
plt.plot(x4, y4, color='b', linewidth=1)
plt.yscale('log')
plt.xlabel('Wavelength ' + '(' + u'\u212B' + ')')
plt.ylabel('log (Flux)')

plt.subplot(112)
#pyplot.savefig('averaged_spectra.eps', format='eps') # saves plot as .eps
plt.savefig('averaged_spectra.pdf', format='PDF')
plt.show()