import numpy as np
import os
import scipy
from astropy.table import Table
import scipy.interpolate as inter

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

for i in range(len(spectra_files)):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
        spectra_name.append(spectra_files[i])
    except ValueError:
        bad_files.append(spectra_files[i])
#        print bad_files

#deredshift data


#scale spectra
wave_min=0  #arbitrary minimum of wavelength range
wave_max=1000000   #arbitrary maximum of wavelength range

for i in range(len(spectra_arrays)):
    spectra = spectra_arrays[i]
#    print spectra
    if (min(spectra[:,0]) > wave_min): #changes minimum wavelength if larger than previous
        wave_min=min(spectra[:,0])
        print spectra_name[i]
    if (max(spectra[:,0]) < wave_max):  #changes maximum wavelength if smaller than previous
        wave_max=max(spectra[:,0])
#        print spectra_name[i]
print wave_min,wave_max


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
	fitted_flux.append(y1)

#sum = fitted_flux[0]
#sums fluxes to be averaged
#for i in range(len(fitted_flux[0])):
#	sum[i] = 0

#for i in range(len(fitted_flux)):
#	sum = sum + fitted_flux[i]


#avg_flux=sum/len(spectra_arrays)	#averages fluxes
#avg_spectrum=Table([wavelength,avg_flux],names=('col1','col2'))	#puts together the average spectrum

#RMS Spectrum

#RMS residual

#plot composite spectrum with RMS Spectrum on top and Residual RMS on bottom
