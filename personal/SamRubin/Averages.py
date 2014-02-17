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
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = supernova()
    if row[2] != None:
        SN.name = row[1]
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




compare_spectrum = SN_Array[0]


"""scale"""
print "Scaling.."
q = 1
for SN in SN_Array[1:]:
    high_wave = 4000
    low_wave = 3000
    if SN.wavelength.min >= low_wave:
        low_wave = SN.wavelength.min
    if SN.wavelength.max <= high_wave:
        high_wave = SN.wavelength.max
    lowindex = np.where(SN.wavelength == low_wave)
    highindex = np.where(SN.wavelength == high_wave)
    print lowindex, "low", highindex, "high"
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
    if SN != SN_Array[0]:
        factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
        scale_factor = np.mean(factors) #"""Maybe make this a weighted average"""
        SN.flux[low:high] *= scale_factor
        SN.error[low:high] *= scale_factor
    print q, "items scaled"
    q += 1

"""composite"""
avg_flux = []
for SN in SN_Array:
    waverange = SN.wavelength[low:high]
    for i in range(len(waverange)):
        average = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
        avg_flux[i] = average
        
        
p1,= plt.plot(waverange, avg_flux)
plt.show()