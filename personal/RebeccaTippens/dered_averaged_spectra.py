from __future__ import division #Future statement to use New-style division
from matplotlib.pyplot as plt
import os #Uses operating system-dependent functionality
import glob #Enables filename pattern matching
import numpy as np 
import scipy as sy 

"""Ask Michael about this: Can't figure out what the initial lo and hi arguments are or exactly how this part is working
I think we should switch to a binary search since we know we're never going to exceed a certain wavelength"""
def find_lo_hi(lo, hi, arr):
    d = arr[:, 0] #Passes the tuple (slice(none, none, none), 0)
    if np.argmin(d) > lo:
        lo = min(d)
    if np.argmax(d) < hi:
        hi = max(d)
    return lo, hi

#Uses the high and low values from find_lo_hi to trim each wavelength array
def trim_data(lo, hi, arr):
    start = np.where(arr == lo)
    end = np.where(arr == hi)
	trimmed = arr[start[0]:end[0]]
	return trimmed

#
path = '../../data/'
files = os.listdir(path)

#Initialize lo and hi to min/max possible
lo = 0.0
hi = float('inf')

for f in files:
    print f
    if f.endswith('.flm'):
        data = np.loadtxt(f)
    else:
        continue
    print data
    lo, hi = find_lo_hi(lo, hi, data)

spectra = {}

for f in files:
    if f.endswith('.flm') or f.endswith('.dat'):
        data = np.loadtxt(f)
    else:
        continue
    trimmed = trim_data(lo, hi, data)
    #scale data by dividing by median
    med = np.median(data[:, 1])
    scaled = np.divide(trimmed, med)
    spectra[f] = trimmed
    print spectra
#-------------- This is the point 
#Change this to average by dividing by the number of files, and generalize the averaging to an actual function
#np.savetxt('dered_averaged_spectra.txt', np.transpose([x3,y3]), fmt="%d %26.18e")
#average flux
avx = g1[:, 0]
#avy = np.divide(np.add(sg1, sg2), 2)
avy = np.divide((sg1+sg2), 2)

#plot trimmed wavelength vs. trimmed/scaled flux
p1 ,= plt.plot(avx, avy, 'c')
p2 ,= plt.plot(g1[:,0], sg1, 'm')
p3 ,= plt.plot(g2[:,0], sg2, 'b')

#Generalize further by making avg an actual fun
def avg(*vals):
	return float(sum(vals))/float(len(vals))

#Plot arrays one by the other
pyplot.plot(x1, y1, color='c', linewidth=1)
pyplot.plot(x2, y2, color='m', linewidth=1)
pyplot.plot(x3, y3, color='#ff6600', linewidth=1)
pyplot.yscale('log')
pyplot.xlabel('Wavelength ' + '(' + u'\u212B' + ')')
pyplot.ylabel('log (Flux)')
pyplot.savefig('dered_averaged_spectra.pdf', format='PDF')
pyplot.show()