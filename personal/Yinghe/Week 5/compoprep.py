# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:20:45 2014

@author: QuantumMonkey
"""
# 
import numpy as np
#import scipy
#import math
import glob
from astropy.table import Table
from astropy.io import ascii
#from astropy.table import Column
import scipy.interpolate as inter
#import sqlite3 as sq3
#import matplotlib.pyplot as plt
import specutils

#Pre-allocate arrays
spectra_files = []
spectra_arrays = []
spectra_name = []
bad_files = []

# read in input spectra

spectra_files=glob.glob("../../../data/cfa/*/*.flm")
                    
#Read in data, store unreadable files (give an error when read in) in a separate array
num0 = 50 # the number of spectra to analyze
for i in range(num0):
    try:
        spectra_arrays.append(Table.read(spectra_files[i],format='ascii'))
        spectra_name.append(spectra_files[i][27:-4])
#        print spectra_name[i]
    except ValueError:
        bad_files.append(spectra_files[i])

num = len(spectra_arrays) # Reset the number to the number of good files
#print num

# dereddening from dust and milky way
# Use model- Fitzpatrick & Massa (2007) extinction model for R_V = 3.1.
#extinction_fm07(wave, a_v=1) 


# deredshift (from old code)
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


# dereddening in supernova rest-frame (from the host galaxy)
# initially assume host galaxy redshift zero





# data interpolation
# pixel size: every 10 As (subject to change)

# still working on putting into one big data structure
wave_min = 1000
wave_max = 20000
wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/10+1)  #creates N equally spaced wavelength values
fitted_flux = []
new = []
#new = Table()
#new['col0'] = Column(wavelength,name = 'wavelength')
for i in range(num):
    new_spectrum=spectra_arrays[i]	#declares new spectrum from list
    new_wave=new_spectrum['col1']	#wavelengths
    new_flux=new_spectrum['col2']	#fluxes
    lower = new_wave[0] # Find the area where interpolation is valid
    upper = new_wave[len(new_wave)-1]
    lines = np.where((new_wave>lower) & (new_wave<upper))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    indata=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    fitted_flux=inter.splev(wavelength,indata)	#fits b-spline over wavelength range
    badlines = np.where((wavelength<lower) | (wavelength>upper))
    fitted_flux[badlines] = 0  # set the bad values to ZERO !!! 
    new = Table([wavelength,fitted_flux],names=('Wavelength','Flux')) # put the interpolated data into the new table
    #newcol = Column(fitted_flux,name = 'Flux')  
    #new.add_column(newcol,index = None)
    
    

    # output data into a file (just for testing, no need to implement)
    output = 'testdata/modified-%s.dat'%(spectra_name[i])
    ascii.write(new, output)
    
    # plot spectra (just for testing, no need to implement)

#    plt.plot(wavelength,fitted_flux)
#    plt.show()

 

