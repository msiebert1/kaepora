import os
import glob
from specutils import extinction as ex
#import astroquery
#from astroquery.ned import Ned
#from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import math
from astropy import units as u
from specutils import Spectrum1D
#import sqlite3 as sq3
#import msgpack


def dered(ebv, wave, flux):
    r_v = 3.1
    flux *= ex.reddening(wave, a_v=ebv*r_v, r_v=3.1, model='f99')
    #need to return just the array!
    return flux


# Data Interpolation

"""
NOTE:
This function inputs three lists containing wavelength, flux and variance.
The output will be a Table with all the fitted values.
You can change the mininum and maximum of wavelength to output,
as well as pixel size in the first few lines.
For the spectra that does not cover the whole range of specified wavelength,
we output the outside values as NAN
"""

from datafidelity import *  # Get inverse variance from the datafidelity outcome

def Interpo (wave, flux, ivar):
    wave_min = 1500
    wave_max = 12000
    dw = 2

    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)
    wavelength = np.arange(math.ceil(wave_min), math.floor(wave_max),
                           dtype=int, step=dw)  # creates N equally spaced wavelength values
    bad_points = []
    inter_flux = []
    inter_ivar = []
    output = []

    lower = wave[0]  # Find the area where interpolation is valid
    upper = wave[-1]

    #ivar = clip(wave, flux, ivar) #clip bad points in flux (if before interpolation)
    ivar = clipmore(wave,flux,ivar)    
    bad_points = clip(wave, flux, ivar)  # if returned bad points range instead of ivar
#    print 'ivar', ivar
#    print 'bad points', bad_points
    #ivar[ivar < 0] = 0 # make sure no negative points

    good_data = np.where((wave >= lower) & (wave <= upper))  #creates an array of wavelength values between minimum and maximum wavelengths from new spectrum

    influx = inter.splrep(wave[good_data], flux[good_data])  # creates b-spline from new spectrum

    inivar = inter.splrep(wave[good_data], ivar[good_data])  # doing the same with the inverse varinces

    inter_flux = inter.splev(wavelength, influx)	 # fits b-spline over wavelength range
    inter_ivar = inter.splev(wavelength, inivar)   # doing the same with errors

#    inter_ivar = clip(wavelength, inter_flux, inter_var) #clip bad points (if do after interpolation)

    # Then the below code (or something similar) would do it (A.S.)
    for wave_tuple in bad_points:
#        print wave_tuple
        zero_points = np.where((wavelength > wave_tuple[0]) & (wavelength < wave_tuple[1]))
        inter_ivar[zero_points] = 0

    inter_ivar[inter_ivar < 0] = 0  # make sure there are no negative points!
    
#    place = np.where((wavelength > 5800.0 ) & (wavelength < 6000.0 ))
#    print inter_ivar[place]

    missing_data = np.where((wavelength < lower) | (wavelength > upper))
    inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
    inter_ivar[missing_data] = float('NaN')

#    print inter_ivar[place]

    output = np.array([wavelength, inter_flux, inter_ivar])  # put the interpolated data into the new table

    return output  # return new table


    # Get the Noise for each spectra ( with input of inverse variance)

def getsnr(flux, ivar):
    sqvar = map(math.sqrt, ivar)
    snr = flux/(np.divide(1.0, sqvar))
    snr_med = np.median(snr)
    return snr_med


def compprep(waves, fluxes, variances, redshift, ebv, dereddened, deredshifted, u_fluxes):
    old_wave = waves
    old_flux = fluxes

    old_wave_r = waves*u.Angstrom 
    old_flux_r = fluxes*u.Unit(u_fluxes)	# fluxes
    old_error = variances  # check if supernovae has error array
    spec1d = Spectrum1D.from_array(old_wave, old_flux)

    old_ivar = genivar(old_wave, old_flux, old_error)  # generate inverse variance
    snr = getsnr(old_flux, old_ivar)
    newdata = []

    if not dereddened:
        new_flux = dered(ebv, old_wave_r, old_flux_r)  # Dereddending (see if sne in extinction files match the SN name)
    if not deredshifted:
        new_wave = old_wave/(1.+redshift)  # Deredshifting

    new_error = old_error  # Placeholder if it needs to be changed
    new_ivar = genivar(new_wave, new_flux, new_error)  # generate new inverse variance
    #var = new_flux*0+1
    newdata = Interpo(new_wave, new_flux, new_ivar)  # Do the interpolation
#    print 'new spectra',newdata
    return newdata, snr
