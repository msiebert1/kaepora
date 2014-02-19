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
compare_spectrum = []

#names = np.loadtxt("", usecols = (0)) #"need to grab all the other data too"
#for row in names:
#    SN = supernova()
#    SN.input = row[0]
#    SN_Array.append(SN)
#print len(names), #"supernovae found"
file_list = []
file_list = glob.glob("../../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../SamRubin/MaxSpectra.dat", dtype = 'str')

# for i in range(len(files)):
# 	file=files[i]
# 	file_list.append(file["col2"])

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading supernovae from database..."
names = []
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = supernova()
    if row[2] != None:
        SN.filename = row[0]
        SN.name = row[1]
        SN.minwave = row[3]
        SN.maxwave = row[4]
        SN_Array.append(SN)
        names.append(SN.name)

print len(SN_Array), "items found"

names,temp = np.unique(names,return_index=True)
SN_Array =  [SN_Array[i] for i in temp]
print len(SN_Array), "unique spectra found"
# rand = np.random.randint(0,len(SN_Array),10)
# all_array = SN_Array
# SN_Array=[]
# for i in rand:
#     tempSN = all_array[i]
#     SN_Array.append(tempSN)
#"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1

for SN in SN_Array:
    k = 1
    for filename in max_light:   
        if SN.filename in filename:
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


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def average(compare_spectrum,SN):
    avg_flux = compare_spectrum.flux
    lowindex = find_nearest(compare_spectrum.wavelength,np.min(SN.wavelength))
    highindex = find_nearest(compare_spectrum.wavelength,np.max(SN.wavelength))
#should be np.where(compare_spectrum.wavelength == np.min(SN.wavelength) if data aligned)
    for i in lowindex+np.arange(highindex-lowindex):
#        avg_flux[i] = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
        try:
            fluxes = np.array([compare_spectrum.flux[i],SN.flux[i]])
            weights = np.array([1./compare_spectrum.error[i], 1./SN.error[i]])    
            avg_flux[i] = np.average(fluxes,weights=weights)
        except IndexError:
            print "Index out of bounds"
            break
    compare_spectrum.flux = avg_flux
# Add residual
    return compare_spectrum

"""scale"""
print "Scaling.."
q = 1
compare_spectrum = SN_Array[0]
for SN in SN_Array:
    low_wave = 4500
    high_wave = 5000
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
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
#    plt.plot(SN.wavelength, SN.flux,label=SN.name)
    print lowindex, "low index", highindex, "high index"
    factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
    scale_factor = np.mean(factors) #"""Maybe make this a weighted average"""
    SN.flux[lowindex:highindex] *= scale_factor
    SN.error[lowindex:highindex] *= scale_factor
    print q, "items scaled"
    compare_spectrum = average(compare_spectrum,SN)
        
    q += 1
"""composite"""

plt.plot(compare_spectrum.wavelength, compare_spectrum.flux,'--',label='Composite')
plt.legend(prop={'size':10})
plt.savefig('Average.png')
plt.show()
