#!/usr/bin/env python -i
'''
Created on Feb 14, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table

class supernova(object):
    """Attributes can be added"""

SN_Array = []

#names = np.loadtxt("", usecols = (0)) #"need to grab all the other data too"
#for row in names:
#    SN = supernova()
#    SN.input = row[0]
#    SN_Array.append(SN)
#print len(names), #"supernovae found"
file_list = []
files = Table.read('../AaronBeaudoin/week4/MaxSpectra.dat',format='ascii')

for i in range(len(files)):
	file=files[i]
	file_list.append(file["col2"])

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading supernovae from database..."
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = supernova()
    if row[2] != None:
        SN.name = row[1]
        SN.minwave = row[3]
        SN.maxwave = row[4]
        SN_Array.append(SN) 
print len(SN_Array), "items found"
        
# rand = np.random.randint(0,len(SN_Array),10)
# all_array = SN_Array
# SN_Array=[]
# for i in rand:
#     tempSN = all_array[i]
#     SN_Array.append(tempSN)
#"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1
for SN in SN_Array[0:50]:
    k = 1
    for filename in file_list:   
        if SN.name in filename:
            data = np.loadtxt(filename)
            SN.wavelength = data[:,0]
            SN.flux = data[:,1]
            error = data[:,2]
            if not all(x == 0.0 for x in error):
                SN.error = error
            break
    print j, 'supernovae fluxed'
    j += 1
print "Done searching. Checking data..."   
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]

print "Done checking. ", len(SN_Array), "items remain"       


def average(compare_spectrum,SN):
    avg_flux = compare_spectrum.flux
    lowindex = np.where(compare_spectrum.wavelength == np.min(SN.wavelength))
    lowindex = lowindex[0]
    highindex = np.where(compare_spectrum.wavelength == np.max(SN.wavelength))
    highindex = highindex[0]

    if len(lowindex) == 0 or len(highindex) == 0: #temperary test
        return compare_spectrum

    for i in lowindex+range(highindex-lowindex):
#        avg_flux[i] = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
        fluxes = np.array([compare_spectrum.flux[i],SN.flux[i]])
        weights = np.array([1./compare_spectrum.error[i], 1./SN.error[i]])
        avg_flux[i] = np.average(fluxes,weights=weights)

    return avg_flux

"""scale"""
print "Scaling.."
q = 1
for SN in SN_Array:
    low_wave = 3500
    high_wave = 4000
#     plt.plot(SN.wavelength, SN.flux,label=SN.name)
    if np.min(SN.wavelength)>= low_wave:
        low_wave = np.min(SN.wavelength)
    if np.max(SN.wavelength) <= high_wave:
        high_wave = np.max(SN.wavelength)
    temp = np.abs(SN.wavelength-low_wave)
    lowindex = np.where(temp == np.min(temp))
    temp = np.abs(SN.wavelength-high_wave)
    highindex = np.where(temp == np.min(temp))
    lowindex = lowindex[0]
    highindex = highindex[0]
    print lowindex, "low", highindex, "high"
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
    plt.plot(SN.wavelength, SN.flux,label=SN.name)

    if SN == SN_Array[0]:
        compare_spectrum = SN
        continue
    factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
    scale_factor = np.mean(factors) #"""Maybe make this a weighted average"""
    SN.flux[lowindex:highindex] *= scale_factor
    SN.error[lowindex:highindex] *= scale_factor
    print q, "items scaled"
    compare_spectrum.flux = average(compare_spectrum,SN)
        
    q += 1
"""composite"""


q = 0

# plt.plot(compare_spectrum.wavelength, compare_spectrum.flux,'--',label='Composite')
plt.legend(prop={'size':10})
plt.show()
