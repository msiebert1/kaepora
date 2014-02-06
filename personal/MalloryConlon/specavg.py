#Written by Celeste, Aaron and Mallory.

import numpy as np
import glob
import scipy
import math
#from astropy.table import Table
import scipy.interpolate as inter
import matplotlib.pyplot as plt

#Pre-allocate arrays
spectra_files = []
spectra_arrays = []
spectra_name = []
bad_files = []

#Start in the data folder, and make a list of each file ending with .flm in each directory


spectra_files=glob.glob("../../data/cfa/*/*.flm")

#Read in data, store unreadable files (give an error when read in) in a separate array
num = 20 # the number of spectra to analyze, change to len(spectra_files) to do all data files in the directory.
for i in range(num):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
        spectra_name.append(spectra_files[i])
    except ValueError:
        bad_files.append(spectra_files[i])


#Deredshift data


#scale/truncate spectra
wave_min=0 #arbitrary minimum of wavelength range
wave_max=1000000 #arbitrary maximum of wavelength range

for i in range(len(spectra_arrays)):
    spectra = spectra_arrays[i]
    if (min(spectra[:,0]) > wave_min): #changes minimum wavelength if larger than previous
        wave_min=min(spectra[:,0])
    # print spectra_name[i]
    if (max(spectra[:,0]) < wave_max): #changes maximum wavelength if smaller than previous
        wave_max=max(spectra[:,0])
# print spectra_name[i]
#print wave_min,wave_max
wavelength = np.linspace(wave_min,wave_max,100) #creates 100 equally spaced wavelength values between the smallest range
fitted_flux=[]

#Generates composite spectrum
for i in range(len(spectra_arrays)):
    new_spectrum=spectra_arrays[i]	#declares new spectrum from list
    new_wave=new_spectrum[:,0]	#wavelengths
    new_flux=new_spectrum[:,1]	#fluxes
    lines=np.where((new_wave>wave_min) & (new_wave<wave_max))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    sm1=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    flux1=inter.splev(wavelength,sm1)	#fits b-spline over wavelength range
    flux1 /= np.median(flux1)
    fitted_flux.append(flux1)
#print fitted_flux

Comp_flux = []
RMS1_flux = []
RMS2_flux = []
for i in range(len(wavelength)):
    fluxa = sum(flux[i] for flux in fitted_flux)
    Comp_flux.append(fluxa/num)
    fluxb = sum((flux[i]-Comp_flux[i])**2 for flux in fitted_flux)
    fluxc = sum((flux[i])**2 for flux in fitted_flux)
    print fluxb,fluxc,Comp_flux[i]
    RMS1_flux.append(Comp_flux[i]+math.sqrt(fluxb/num))
    RMS2_flux.append(Comp_flux[i]-math.sqrt(fluxb/num))
#print RMS_flux

#plot composite spectrum
plot1 = plt.plot(wavelength,Comp_flux,label = 'comp')
plot2 = plt.plot(wavelength,RMS1_flux,label = 'rms+')
plot3 = plt.plot(wavelength,RMS2_flux,label = 'rms-')
legend = plt.legend(loc='upper right', shadow=True)
plt.xlim(wave_min,wave_max)
plt.xlabel('Observed Wavelength ($\AA$)') #Change to Rest Wavelength when redshift is figured out
plt.ylabel('Scaled Flux')
#plt.yscale('log')
plt.show()
plt.savefig('compositerms.png')

#RMS Spectrum

#RMS residual

#plot composite spectrum with RMS Spectrum on top and Residual RMS on bottom
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
ax0.plot()
ax0.set_title('RMS')
ax0.set_ylabel('')
ax1.plot()
ax1.set_title('Residual RMS')
ax1.set_yscale('log')
ax1.set_ylabel('')
plt.subplots_adjust(hspace=0.3)
plt.show()
#plt.savefig('')