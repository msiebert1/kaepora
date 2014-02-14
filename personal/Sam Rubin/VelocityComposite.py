'''
Created on Feb 9, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sqlite3 as sq3
import scipy
from scipy import interpolate as intp
import math


#creates a new class for storing all the relevant data
class Supernovae(object):
    """HELLO"""
    
#connect to the database
con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()

#Creates an array containing only the used information, which can be easily referenced 
SN_Array = []
print "Reading supernovae data..."
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = Supernovae()
    SN.filename = row[0]
    SN.name = row[1]
    if row[2] == None:
        pass
    else:
        SN.redshift = float(row[2])
    SN.minwave = row[3]
    SN.maxwave = row[4]
    SN_Array.append(SN)
print len(SN_Array)

#open the velocity file, make it useable, add the needed data to the array
v_data = np.loadtxt('foley_master_data', dtype={'names': ('name', 'redshift', 'v', 'dv'), 'formats': ('S8', 'f', 'f', 'f')}, skiprows = 1, usecols = (0,1,2,3))
print "Reading velocities..."
for row in v_data:
    for i in range(len(SN_Array)):
        if row[0] == SN_Array[i].name:
            SN_Array[i].redshift = float(row[1])
            SN_Array[i].v_si = row[2]
            SN_Array[i].dv_si = row[3]
print len(v_data)
#Takes out files that didn't have velocity data
print "Cleaning..."       
for SN in SN_Array:
    if not (hasattr(SN, 'v_si')):
        SN_Array.remove(SN)
print len(SN_Array)
print "Globbing files"
spectra_files = glob.glob("../../data/cfa/*/*.flm")
#Adds flux to the stuff
print "Grabbing flux data..."
SN_Array = SN_Array[0:9] #makes this way more manageable
j = 1
for SN in SN_Array:
    print j
    j += 1
    for filename in spectra_files[0:100]:
        try:
            data = np.loadtxt(filename)
            SN.wavelength = data[:,0]
            SN.flux = data[:,1]
            SN.flux_error = data[:,2]
        except ValueError:
            print "Invalid file!"
print len(SN_Array)

#wavelength linspace
wave_min = math.floor(max(SN.minwave for SN in SN_Array))
wave_max = math.floor(min(SN.maxwave for SN in SN_Array))
waverange = np.linspace(wave_min, wave_max, (wave_max-wave_min)*10)
waverange = np.array(waverange)

print "De-redshifting..."
#De-redshift, scale, and interpolate
SN_Array_high = []
SN_Array_low = []
for i in range(len(SN_Array)):
    SN_Array[i].wavelength /= (1 + SN_Array[i].redshift)
    spline_rep = intp.splrep(SN_Array[i].wavelength, SN_Array[i].flux)
    new_flux = intp.splev(waverange, spline_rep)
    new_flux /= np.median(new_flux)
    SN_Array[i].wavelength = waverange
    SN_Array[i].flux = new_flux
    if SN_Array[i].v_si <= -11:
        SN_Array_high.append(SN_Array[i])
    else:
        SN_Array_low.append(SN_Array[i])
    
#generate composite spectra 
comp_flux_high = []
comp_flux_low = []

#Weighted average
print "Averaging..."  
for i in range(waverange):
    weights = np.array(1/(SN.flux_error[i]**2) for SN in SN_Array_high)
    fluxes = np.array(SN.flux[i] for SN in SN_Array_high)
    comp_flux_high.append(sum(weights*fluxes)/sum(weights))
for i in range(waverange):
    weights = np.array(1/(SN.flux_error[i]**2) for SN in SN_Array_low)
    fluxes = np.array(SN.flux[i] for SN in SN_Array_low)
    comp_flux_low.append(sum(weights*fluxes)/sum(weights))
print "Done!"

#Plot that shit
p1,=plt.plot(waverange, comp_flux_high)
p2,=plt.plot(waverange, comp_flux_low)
plt.xlabel('Wavelength [A]')
plt.ylabel('Scaled Flux')
plt.yscale('log')
plt.legend([p1,p2],['High V_Si', 'Low V_Si'],2)
plt.savefig('Composite_by_Velocity.png')
plt.show()

con.close()