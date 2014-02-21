'''
Created on Feb 9, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math


"""creates a new class for storing all the relevant data"""
class Supernovae(object):
    """Empty now, will be filled"""

"""Creates an array containing only the used information, which can be easily referenced """
SN_Array = []
"""open the velocity file, make it useable, add the needed data to the array"""
print "Reading velocities..."
v_data = np.loadtxt('foley_master_data', dtype={'names': ('name', 'redshift', 'v', 'dv'), 'formats': ('S8', 'f8', 'f8', 'f8')}, skiprows = 1, usecols = (0,1,2,3))
for row in v_data[0:100]:
    SN = Supernovae()
    SN.name = row[0]
    SN.redshift = row[1]
    SN.v_si = row[2]
    SN.dv_si = row[3]
    SN_Array.append(SN)
print len(v_data), "velocities found"

#read in data structure
"""connect to the database"""
con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
"""Now deal with the database"""
print "Reading supernovae from database..."
for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, MinWave, MaxWave FROM Supernovae'):
        if row[1] == SN.name:
            SN.filename = row[0]
            SN.minwave = row[2]
            SN.maxwave = row[3]
            break    
print len(SN_Array), "items found"
"""Takes out files that didn't have velocity data"""
print "Cleaning..."       
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'filename')]
print len(SN_Array), "items remain"

"""Grab list of every spectrum"""
print "Globbing file list"
spectra_files = glob.glob('../../data/cfa/*/*.flm')

"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1
for SN in SN_Array:
    k = 1
    for filename in spectra_files:   
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

"""wavelength linspace"""
wave_min = math.floor(max(SN.minwave for SN in SN_Array))
wave_max = math.floor(min(SN.maxwave for SN in SN_Array))
waverange = np.linspace(wave_min, wave_max, (wave_max-wave_min)*10)
waverange = np.array(waverange)

print "De-redshifting, interpolating, and scaling..."
"""De-redshift, scale, and interpolate"""
SN_Array_high = []
SN_Array_low = []
high_v = -12
for SN in SN_Array:
    SN.wavelength = np.divide(SN.wavelength, 1 + SN.redshift)
    try:
        spline1 = intp.splrep(SN.wavelength, SN.flux)
        spline2 = intp.splrep(SN.wavelength, SN.error)
    except ValueError:
        print "Invalid data found"
    new_flux = intp.splev(waverange, spline1)
    new_flux /= np.median(new_flux)
    new_error = intp.splev(waverange, spline2)
    SN.wavelength = waverange
    SN.flux = new_flux
    SN.error = new_error
    if SN.v_si <= high_v:
        SN_Array_high.append(SN)
    else:
        SN_Array_low.append(SN)
print "Done with that. ", len(SN_Array_high), "high-V supernovae found. ", len(SN_Array_low), "low-V supernovae found."
    
"""generate composite spectra by weighted average"""
print "Averaging..."  
comp_flux_high = []
comp_flux_low = []
for i in range(len(waverange)):
    try:
        flux = np.sum(SN.flux[i]/(SN.error[i])**2 for SN in SN_Array_high if SN.error[i] != 0)/np.sum(1/(SN.error[i])**2 for SN in SN_Array_high if SN.error[i] != 0 and SN.error[i] != None)
    except ValueError:
        print "Invalid data found"
    comp_flux_high.append(flux)
print "High-V supernovae have been averaged."
for i in range(len(waverange)):
    try:
        flux = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array_low if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array_low if SN.error[i] != 0 and SN.error[i] != None)
    except ValueError:
        print "Invalid data found"        
    comp_flux_low.append(flux)
print "Low-V supernovae have been averaged, too."

"""Plot that shit"""
print "Plotting..."
p1,=plt.plot(waverange, comp_flux_high)
p2,=plt.plot(waverange, comp_flux_low)
plt.xlabel('Wavelength [A]')
plt.ylabel('Scaled Flux')
"""plt.yscale('log')"""
plt.legend([p1,p2],['High V_Si', 'Low V_Si'],2)
plt.savefig('Composite_by_Velocity.png')
plt.show()
print "Done!"

con.close()