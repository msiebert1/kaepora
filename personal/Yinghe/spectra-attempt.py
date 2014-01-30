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
# not sure the input form of the redshift of each supernova. Could also use an array.
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

# Function: import into an array, and averaging
def ImpArr(d,z) :
    walen = []
    flux = []
    fluxavg = []            
    for i in range(len(d)) :
        if d[i][1] > 0 :
            walen.append(d[i][0]/(1+z))
            flux.append(log(d[i][1]))
#    print flux
    avg = float(sum(flux))/len(flux)
    print avg
    for i in range(len(flux)) :
        fluxavg.append(flux[i]-avg)
#        print fluxavg[i]
    return [walen,flux,fluxavg]  

# main function (can be updated)        
data1 = 'sn2011by-hst+lick_sdx.txt'
data2 = 'sn2011fe-visit3-hst_sdx.txt'
#for n in range(2):
d1 = ReadIn(data1)
arr1 = ImpArr(d1,z1)
d2 = ReadIn(data2)
arr2 = ImpArr(d2,z2)
#     
mini = max(arr1[0][0],arr2[0][0])
maxi = min(arr1[0][len(arr1[0])-1],arr2[0][len(arr2[0])-1]) 
#print mini,maxi

# plotting spectra
f1 = plt.figure(1)
plot1 = plt.plot(arr1[0],arr1[1])
plot2 = plt.plot(arr2[0],arr2[1])
plt.xlabel('Rest Wavelength [A]')
plt.ylabel('Flux')
plt.xlim(mini,maxi)
#plt.ylim()
f1.show()
f1.savefig('test1.png')

# plotting average spectra
f2 = plt.figure(2)
plot3 = plt.plot(arr1[0],arr1[2])
plot4 = plt.plot(arr2[0],arr2[2])
l = plt.axhline(y=0,color = 'r')
plt.xlabel('Rest Wavelength')
plt.xlim(mini,maxi)
plt.ylabel('Average Flux')
f2.show()
f2.savefig('test2.png')
