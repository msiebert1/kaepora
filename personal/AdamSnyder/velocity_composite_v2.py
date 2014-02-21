import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
import math
from scipy import interpolate
import sqlite3 as sq3

# connect to the database

# Define a Supernovae and intialize list
class Supernovae(object):
    """attributes: Filename, name, redshift, flux-wavelength array, MinWave, MaxWave"""    

SN_Array = [] 
SN_High = []
SN_Low = []
high_v_Si = -11

#Load supernovae to analyze from data file 

v_data = np.loadtxt('foley_master_data', dtype={'names': ('name', 'z', 'v', 'dv'), 'formats': ('S8', 'f', 'f', 'f')}, skiprows = 1, usecols = (0, 1, 2, 3))

print "Reading data file..."

for line in v_data:
    SN = Supernovae()
    SN.name = line[0]
    if line[1] >=0.0:
        SN.redshift = line[1] #check for anomalous redshift
    SN.v_Si = line[2]
    SN.dv_Si = line[3]
    SN_Array.append(SN)

print "Done reading.", len(SN_Array), "Supernovae read. Opening SQL database..."

# Connect to SQL database

con = sq3.connect('../MichaelSchubert/Sne.db')
cur = con.cursor()

for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, MinWave, MaxWave FROM Supernovae'):
        if row[1] == SN.name:
            SN.filename = row[0]
            SN.MinWave = row[2]
            SN.MaxWave = row[3]
            break

print "Cleaning array..."

SN_Array = [SN for SN in SN_Array if hasattr(SN, 'filename') and hasattr(SN, 'redshift')]

print len(SN_Array), "supernovae read."
print "Reading lists of spectra files to find the fluxes"

spectra_files = glob.glob("../../data/cfa/*/*.flm")

print "Searching for matching flux data."

for SN in SN_Array:
    for filename in spectra_files:
        if SN.filename in filename:
            data = np.loadtxt(filename)
            SN.wavelength = data[:, 0]
            SN.flux = data[:, 1]
            error = data[:, 2]
            if not all(x == 0.0 for x in error):
                SN.error = error
            break

SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]

print "Done searching.", len(SN_Array), "supernovae read."

wave_min = math.floor(max(SN.MinWave for SN in SN_Array))
wave_max = math.floor(min(SN.MaxWave for SN in SN_Array))
wavelength = scipy.linspace(wave_min, wave_max, (wave_max-wave_min)*10)

print "Formatting data."

for SN in SN_Array:
    SN.wavelength = np.divide(SN.wavelength, 1 + SN.redshift)
    spline_rep1 = interpolate.splrep(SN.wavelength, SN.flux)
    spline_rep2 = interpolate.splrep(SN.wavelength, SN.error)
    new_flux = interpolate.splev(wavelength, spline_rep1)
    new_flux /= np.median(new_flux)
    new_error = interpolate.splev(wavelength, spline_rep2)
    SN.wavelength = wavelength
    SN.flux = new_flux
    SN.error = new_error
    if SN.v_Si <= high_v_Si:
        SN_High.append(SN)
    else:
        SN_Low.append(SN)

print "Done formatting data."
print len(SN_High), "high velocity supernovae"
print len(SN_Low), "low velocity supernovae"

Comp_flux_high = []
Comp_flux_low = []

for i in range(len(wavelength)):
    flux = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_High if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_High if SN.error[i] != 0)
    Comp_flux_high.append(flux)

for i in range(len(wavelength)):
    flux = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Low if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Low if SN.error[i] != 0)
    Comp_flux_low.append(flux)

plt.figure(1)
plt.subplot(211)
plt.plot(wavelength, Comp_flux_high, 'b', wavelength, Comp_flux_low, 'r')
plt.title('High/Low Velocity Composite Fluxes')
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux')
plt.subplot(212)
#plt.plot(stuff)
plt.title('Residual')
plt.xlabel('Wavelength (A)')
plt.ylabel('Percantage')
plt.show()
