# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:54:54 2014

@author: QuantumMonkey
"""
import numpy as np
import scipy.interpolate as inter
from math import floor, ceil
from astropy.table import Table



# the interpolation function

def Interpo(spectra) :
    wave_min = 1000
    wave_max = 20000
    pix = 2
    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)  #creates N equally spaced wavelength values
    wavelength = np.arange(ceil(wave_min), floor(wave_max), dtype=int, step=pix)
    fitted_flux = []
    fitted_error = []
    new = []
    #new = Table()
    #new['col0'] = Column(wavelength,name = 'wavelength')
    new_spectrum=spectra	#declares new spectrum from list
    new_wave=new_spectrum[:,0]	#wavelengths
    new_flux=new_spectrum[:,1]	#fluxes
    new_error=new_spectrum[:,2]   #errors
    lower = new_wave[0] # Find the area where interpolation is valid
    upper = new_wave[len(new_wave)-1]
    lines = np.where((new_wave>lower) & (new_wave<upper))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    indata=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
    inerror=inter.splrep(new_wave[lines],new_error[lines]) # doing the same with the errors
    fitted_flux=inter.splev(wavelength,indata)	#fits b-spline over wavelength range
    fitted_error=inter.splev(wavelength,inerror)   # doing the same with errors
    badlines = np.where((wavelength<lower) | (wavelength>upper))
    fitted_flux[badlines] = 0  # set the bad values to ZERO !!! 
    new = Table([wavelength,fitted_flux],names=('col1','col2')) # put the interpolated data into the new table    
    #newcol = Column(fitted_flux,name = 'Flux')  
    #new.add_column(newcol,index = None)
    return new