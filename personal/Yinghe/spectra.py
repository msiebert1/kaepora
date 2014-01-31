# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:56:57 2014

@author: QuantumMonkey
"""
#!/usr/bin/env/python

import math
import numpy as np
import matplotlib.pyplot as plt

#redshifts here
# not sure how to use
z1 = 0.002843
z2 = 0.000804

# Function: reading data

def ReadIn(file):
    f = open(file,'r')
    d = []
    for line in f.readlines():
        line=line.decode('utf-8-sig')
        line=line.encode('UTF-8')
    #    print line.split()
        d.append([float(value) for value in line.split()])    
    f.close()
    return d

# Function: import into an array, and find the mean value
def ImpArr(d) :
    walen = []
    flux = []           
    for i in range(len(d)) :
        walen.append(d[i][0])
        flux.append(d[i][1])
    med = np.median(flux)
    print med
    flux = flux/med
#    print flux
#    avg = float(sum(flux))/len(flux)
#    print avg
#    for i in range(len(flux)) :
#        fluxavg.append(flux[i]-avg)
#        print fluxavg[i]
    return [walen,flux]
    

# main function (can be updated)        
data1 = 'sn2011by-hst+lick_sdx.txt'
data2 = 'sn2011fe-visit3-hst_sdx.txt'
d1 = ReadIn(data1)
arr1 = ImpArr(d1)
d2 = ReadIn(data2)
arr2 = ImpArr(d2)
#  cut and average the spectra   
# For two-case only. Still working on how to do multi-cases...
mini = max(arr1[0][0],arr2[0][0])
maxi = min(arr1[0][len(arr1[0])-1],arr2[0][len(arr2[0])-1]) 
#print mini,maxi
for i in range(len(arr1[0])) : # what if the range is overlapped?
    if arr1[0][i] == mini : 
        minum = i
    if arr1[0][i] == maxi : 
        manum = i
#print minum,manum
avwalen = []
avflux = []
avwalen= arr1[0][minum:manum+1]
#print avwalen
for k in range(len(avwalen)) :
    avflux.append((arr1[1][k]+arr2[1][k])/2)
#print avfluxnorm

# plotting spectra
f1 = plt.figure(1)
ax = f1.add_subplot(111)
plot1 = ax.plot(arr1[0],arr1[1], label = 'SN2011BY')
plot2 = ax.plot(arr2[0],arr2[1], label = 'SN2011FE')
plot3 = ax.plot(avwalen,avflux, label = 'average')
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux')
plt.xlim(mini,maxi)
ax.set_yscale('log')
#plt.ylim()
legend = ax.legend(loc='lower right', shadow=True)
f1.show()
f1.savefig('spectra.png')

# plotting normalized spectra
"""f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
plot4 = ax2.plot(arr1[0],arr1[2],label = 'SN2011BY')
plot5 = ax2.plot(arr2[0],arr2[2], label = 'SN2011FE')
plot6 = ax2.plot(avwalen,avfluxnorm, label = 'average')
#l = plt.axhline(y=0,color = 'r')
plt.xlabel('Wavelength [A]')
plt.xlim(mini,maxi)
ax2.set_yscale('log')
plt.ylabel('Normalized Flux')
f2.show()
f2.savefig('spectra2.png')"""