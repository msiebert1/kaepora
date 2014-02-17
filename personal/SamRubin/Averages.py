'''
Created on Feb 14, 2014

@author: Alien
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math

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
file_list = glob.glob("../../data/cfa/*/*.flm")

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading supernovae from database..."
for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
        SN = supernova()
        if row[2] != None:
            SN.filename = row[0]
            SN.minwave = row[3]
            SN.maxwave = row[4]
            SN_Array.append(SN) 
print len(SN_Array), "items found"
        
#"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1
for SN in SN_Array[0:100]:
    k = 1
    for filename in file_list:   
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


high_wave = 4000
low_wave = 3000

compare_spectrum = SN_Array[0]


"""scale"""
print "scaling.."
for SN in SN_Array[1:]:
    low = SN.wavelength.index(3000)
    high = SN.wavelength.index(4000)
    if SN.wavelength.min >= low_wave:
        low_wave = SN.wavelength.min
        low = SN.wavelength.index(low_wave)
    if SN.wavelength.max <= high_wave:
        high_wave = SN.wavelength.max
        high = SN.wavelength.index(high_wave)
    SN.flux /= np.median(SN.flux)
    if SN != SN_Array[0]:
        factors = compare_spectrum.flux[low:high] / SN.flux[low:high]
        scale_factor = np.mean(factors) #"""Maybe make this a weighted average"""
        SN.flux[low:high] *= scale_factor
        SN.error[low:high] *= scale_factor

"""composite"""
avg_flux = []
for SN in SN_Array:
    waverange = SN.wavelength[low:high]
    for i in range(len(waverange)):
        average = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
        avg_flux[i] = average
        
        
p1,= plt.plot(waverange, avg_flux)
plt.show()