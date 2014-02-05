
#Measure the scatter for the composite spectrum.  
#Make a two-panel plot which shows the composite spectrum with the RMS spectrum around it on top and the residual RMS spectrum on bottom.

import numpy as np
import os
import scipy
#from astropy.table import Table
import scipy.interpolate as inter
import matplotlib.pyplot as plt

#Pre-allocate arrays
spectra_files = []
spectra_arrays = []
spectra_name = []
bad_files = []

#Start in the data folder, and make a list of each file ending with .flm in each directory

path = '../../data/cfa/'

for dirs,subdirs,files in os.walk(path):
    for file in files:
        if file.endswith('.flm'):
            spectra_files.append(os.path.join(dirs,file))

#Read in data, store unreadable files (give an error when read in) in a separate array

for i in range(100):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
        spectra_name.append(spectra_files[i])
    except ValueError:
        bad_files.append(spectra_files[i])
#        print bad_files

#deredshift data


#scale/truncate spectra
wave_min=0  #arbitrary minimum of wavelength range
wave_max=1000000   #arbitrary maximum of wavelength range

for i in range(len(spectra_arrays)):
    spectra = spectra_arrays[i]
    if (min(spectra[:,0]) > wave_min): #changes minimum wavelength if larger than previous
        wave_min=min(spectra[:,0])
#        print spectra_name[i]
    if (max(spectra[:,0]) < wave_max):  #changes maximum wavelength if smaller than previous
        wave_max=max(spectra[:,0])
#        print spectra_name[i]
#print wave_min,wave_max
wavelength = np.linspace(wave_min,wave_max,100)  #creates 100 equally spaced wavelength values between the smallest range
fitted_flux=[]
#generates composite spectrum
for i in range(len(spectra_arrays)):
    new_spectrum=spectra_arrays[i]	#declares new spectrum from list
    new_wave=new_spectrum[:,0]	#wavelengths
    new_flux=new_spectrum[:,1]	#fluxes
    lines=np.where((new_wave>wave_min) & (new_wave<wave_max))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    sm1=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    y1=inter.splev(wavelength,sm1)	#fits b-spline over wavelength range
    y1 /= np.median(y1)
    fitted_flux.append(y1)
    
Comp_flux = []   
for i in range(len(wavelength)): # For each wavelength, sum associated flux for each spectra and average 
    fluxa = sum(flux[i] for flux in fitted_flux)
    Comp_flux.append(fluxa/float(len(wavelength)))     
#print Comp_flux
    
#plot composite spectrum
plt.plot(wavelength,Comp_flux)
plt.xlim(wave_min,wave_max)
plt.xlabel('Wavelength')
plt.ylabel('Flux')
#plt.yscale('log')
plt.show()
plt.savefig('composite.png')

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
#plt.show()
#plt.savefig('')
