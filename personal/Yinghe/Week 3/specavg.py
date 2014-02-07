
import numpy as np
import os
import scipy
import math
import glob
from astropy.table import Table
import scipy.interpolate as inter
import matplotlib.pyplot as plt

#Pre-allocate arrays
spectra_files = []
spectra_arrays = []
spectra_name = []
bad_files = []

#Start in the data folder, and make a list of each file ending with .flm in each directory

path = '../../../data/cfa/'

"""for dirs,subdirs,files in os.walk(path):
    for file in files:
        if file.endswith('.flm'):
            spectra_files.append(os.path.join(dirs,file))"""
            
#Faster ways            
spectra_files=glob.glob("../../../data/cfa/*/*.flm")
                    
#Read in data, store unreadable files (give an error when read in) in a separate array
num0 = 50 # the number of spectra to analyze
for i in range(num0):
    try:
        spectra_arrays.append(Table.read(spectra_files[i],format='ascii'))
        spectra_name.append(spectra_files[i])
    except ValueError:
        bad_files.append(spectra_files[i])

num = len(spectra_arrays) # Reset the number to the number of good files

#deredshift data
parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
sn_name = parameters['col1']
sn_z = parameters['col2']
#old_spectrum = []
for i in range(num):
    old_spectrum = spectra_arrays[i]
    z=0
    file_name = spectra_name[i]
    for j in range(len(sn_name)):
        if sn_name[j] in file_name:
		z=sn_z[j]
    lambda_obs =old_spectrum['col1']
    lambda_emit= lambda_obs/(1+z)
    spectra_arrays[i]=Table([lambda_emit,old_spectrum['col2']],names=('col1','col2'))

#scale/truncate spectra
wave_min=0  #arbitrary minimum of wavelength range
wave_max=1000000   #arbitrary maximum of wavelength range

for i in range(len(spectra_arrays)):
    spectra = spectra_arrays[i]
    if (min(spectra['col1']) > wave_min): #changes minimum wavelength if larger than previous
        wave_min=min(spectra['col1'])
#        print spectra_name[i]
    if (max(spectra['col1']) < wave_max):  #changes maximum wavelength if smaller than previous
        wave_max=max(spectra['col1'])
#        print spectra_name[i]
#print wave_min,wave_max

wavelength = np.linspace(wave_min,wave_max,100)  #creates 100 equally spaced wavelength values between the smallest range
fitted_flux=[]

#generates composite spectrum
for i in range(num):
    new_spectrum=spectra_arrays[i]	#declares new spectrum from list
    new_wave=new_spectrum['col1']	#wavelengths
    new_flux=new_spectrum['col2']	#fluxes
    lines=np.where((new_wave>wave_min) & (new_wave<wave_max))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    sm1=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    y1=inter.splev(wavelength,sm1)	#fits b-spline over wavelength range
    y1 /= np.median(y1)
    fitted_flux.append(y1)
#print fitted_flux

# Calculating RMS and Residuals    
Comp_flux = []   
RMS1_flux = []
RMS2_flux = []
scatter = []
for i in range(len(wavelength)): 
    fluxa = sum(flux[i] for flux in fitted_flux)
    Comp_flux.append(fluxa/num)
    fluxb = sum((flux[i]-Comp_flux[i])**2 for flux in fitted_flux)
    RMS1_flux.append(Comp_flux[i]+math.sqrt(fluxb/num))
    RMS2_flux.append(Comp_flux[i]-math.sqrt(fluxb/num))
    scatter.append(math.sqrt(fluxb/num)/Comp_flux[i])
#    if RMS1_flux[i] > 10000 : print ERROR
#print Comp_flux,RMS1_flux
    

#plot composite spectrum with RMS Spectrum on top and Residual RMS on bottom
fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
plot1 = ax0.plot(wavelength,Comp_flux,label = 'comp')
plot2 = ax0.plot(wavelength,RMS1_flux,label = 'rms+')
plot3 = ax0.plot(wavelength,RMS2_flux,label = 'rms-')
legend = ax0.legend(loc='upper right', shadow=True)
ax0.set_xlim(wave_min,wave_max)
ax0.set_title('RMS spectrum')
ax0.set_ylabel('Flux')
ax1.plot(wavelength,scatter,'o-')
ax1.set_ylabel('Residual')
ax1.set_ylim(0,1)
legend = ax1.legend(loc='upper right', shadow=True)
#ax1.set_yscale('log')
#ax1.set_ylabel('')
plt.subplots_adjust(hspace=0.1)
plt.show()
plt.savefig('all.png')
