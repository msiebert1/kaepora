# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:20:45 2014

@author: QuantumMonkey
"""
# 
import numpy as np
import os
import scipy
import math
import glob
from astropy.table import Table
import scipy.interpolate as inter
import matplotlib.pyplot as plt


# read in input spectra
#Pre-allocate arrays
spectra_files = []
spectra_arrays = []
spectra_name = []
bad_files = []

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

# dereddening from dust and milky way - use Patrick 




# deredshift
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





# interpolation
# pixel - size 
wave_min = 1000
wave_max = 20000
wavelength = np.linspace(wave_min,wave_max,10000)  #creates 100 equally spaced wavelength values between the smallest range
fitted_flux=[]

for i in range(num):
    new_spectrum=spectra_arrays[i]	#declares new spectrum from list
    new_wave=new_spectrum['col1']	#wavelengths
    new_flux=new_spectrum['col2']	#fluxes
    lines=np.where((new_wave>wave_min) & (new_wave<wave_max))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    sm1=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    y1=inter.splev(wavelength,sm1)	#fits b-spline over wavelength range
    fitted_flux.append(y1)

# varaiance- inverse square of error 

