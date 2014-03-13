# -*- coding: utf-8 -*-
"""
Created on Thu Mar 06 23:14:23 2014

@author: QuantumMonkey
"""

import glob
import numpy as np
from math import floor,ceil

#list of files
spectra_files = glob.glob ('../../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []

junk_data = []

new_spectra = []

#number of spectra to modify
num = 2

#get data, pathnames
for i in range(num):
	try:
        	spectra_data.append(np.loadtxt(spectra_files[i]))
        	file_path.append(spectra_files[i][14:-4])
		#print file_path
             
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)

#import rebinning
from rebinning import *

for i in range(num) :  
    wave = spectra_data[i][:,0]	#wavelengths
    flux = spectra_data[i][:,1]	#fluxes
    wave_min = 1000
    wave_max = 20000
    pix = 2
#    nwave = np.arange(ceil(wave_min), floor(wave_max), dtype=int, step=pix) #creates N equally spaced wavelength values
    nwave = spectra_data[i+1][:,0]
    nflux = spectra_data[i+1][:,1]
#    print nwave,nflux
    new_spectra.append(womashrebin(wave, flux, nwave, nflux))
    print new_spectra
    
    
    
