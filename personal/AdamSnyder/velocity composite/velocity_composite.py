import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
import math
from scipy import interpolate
import sqlite3 as sq3

# connect to the database

# Define a Supernovae 
class Supernovae(object):
    """attributes: Name, redshift, flux-wavelength array, MinWave, MaxWave"""    

con = sq3.connect('../MichaelSchubert/Sne.db')
cur = con.cursor()
high_v = -11.0

SN_Array = []

for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = Supernovae()
    SN.filename = row[0]
    SN.name = row[1]
    SN.redshift = row[2]
    SN.MinWave = row[3]
    SN.MaxWave = row[4]
    SN_Array.append(SN)

print len(SN_Array)

v_data = np.loadtxt('foley_master_data', dtype={'names': ('name', 'z', 'v', 'dv'), 'formats': ('S8', 'f', 'f', 'f')}, skiprows = 1, usecols = (0, 1, 2, 3))


for data in v_data:
    for i in range(len(SN_Array)):
        if data[0] == SN_Array[i].name:
            SN_Array[i].vel = data[2]
            SN_Array[i].dv = data[3]

SN_Array = [SN for SN in SN_Array if hasattr(SN, 'vel')]

spectra_files = glob.glob("../../data/cfa/*/*.flm")

j = 0

for SN in SN_Array[0:10]:
    j += 1
    for filename in spectra_files:
        if SN.filename in filename:
            data = np.loadtxt(filename)
            SN.wavelength = np.array(data[:, 0])
            SN.flux = np.array(data[:, 1])
            SN.error = np.array(data[:, 2])
            
            print 'You have finished {0}%\r',format(j/100.0),

SN_Array = [SN for SN in SN_Array if hasattr(SN, 'flux')]
print len(SN_Array)

SN_High = []
SN_Low = []

wave_min = math.floor(max(SN.MinWave for SN in SN_Array))
wave_max = math.floor(min(SN.MaxWave for SN in SN_Array))
wavelength = scipy.linspace(wave_min, wave_max, (wave_max-wave_min)*10)

for i in range(len(SN_Array)):
    SN_Array[i].wavelength = np.divide(SN_Array[i].wavelength, 1 + SN_Array[i].redshift)
    spline_rep = interpolate.splrep(SN_Array[i].wavelength, SN_Array[i].flux)
    new_flux = interpolate.splev(wavelength, spline_rep)
    new_flux /= np.median(new_flux)
    SN_Array[i].wavelength = wavelength
    SN_Array[i].flux = new_flux
    if SN_Array[i].vel <= high_v:
       SN_High.append(SN_Array[i])
    else:
       SN_Low.append(SN_Array[i])

# Generate composite spectra by averaging

Comp_flux_high = []
Comp_flux_low = []

for i in range(len(wavelength)):
    weights = np.array(list(1/(supernova.error[i]**2) for supernova in SN_High))
    fluxes = np.array(list(supernova.flux[i] for supernova in SN_High))
    Comp_flux_high.append(sum(weights*fluxes)/sum(weights))

for i in range(len(wavelength)):
    weights = np.array(1/(supernova.error[i]**2) for supernova in SN_Low)
    fluxes = np.array(supernova.flux[i] for supernova in SN_Low)
    Comp_flux_low.append(sum(weights*fluxes)/sum(weights))

plt.plot(wavelength, Comp_flux_high, 'b', wavelength, Comp_flux_low, 'r')
plt.title('High/Low Velocity Composite Fluxes')
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux')
plt.show()

"""
res_flux = []

for i in range(num_spectra):
    res_flux.append(np.array(new_spectra[i][1]-Comp_spectra[1, :]))

res_flux = np.array(res_flux)
RMS = np.sqrt(np.mean(res_flux*res_flux, axis = 0))
flux_err_pos = Comp_spectra[1, :] + RMS
flux_err_neg = Comp_spectra[1, :] - RMS
scatter = RMS / Comp_spectra[1, :]

# Plot composite spectra (Jut need to figure out what to plot, and look up to make single figures)

plt.figure(1)
plt.subplot(211)
plt.plot(Comp_spectra[0, :], Comp_spectra[1, :], 'b' , Comp_spectra[0, :], flux_err_pos, 'r', Comp_spectra[0, :], flux_err_neg, 'g')
plt.title('Composite Spectra and (+/-) RMS')
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux')
plt.subplot(212)
plt.plot(wavelength, scatter)
plt.title('Residual')
plt.xlabel('Wavelength (A)')
plt.ylabel('Percantage')
plt.show()

"""
