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
from datafidelity import *  # Get variance from the datafidelity outcome

def Interpo (wave, flux, variance) :
    wave_min = 1500
    wave_max = 12000
    dw = 2

    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)
    wavelength = np.arange(math.ceil(wave_min), math.floor(wave_max), dtype=int, step=dw) #creates N equally spaced wavelength values
    inter_flux = []
    inter_var  = []
    output     = []

    lower = wave[0] # Find the area where interpolation is valid
    upper = wave[-1]

    good_data = np.where((wave >= lower) & (wave <= upper))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum

    influx = inter.splrep(wave[good_data], flux[good_data])	#creates b-spline from new spectrum

    invar  = inter.splrep(wave[good_data], variance[good_data]) # doing the same with the errors

    inter_flux = inter.splev(wavelength, influx)	#fits b-spline over wavelength range
    inter_var  = inter.splev(wavelength, invar)   # doing the same with errors

    new_inter_var = clip(wavelength, inter_flux, inter_var)
    new_inter_var[new_inter_var < 0] = 0

    missing_data = np.where((wavelength < lower) | (wavelength > upper))
    inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
    new_inter_var[missing_data] =  float('NaN')

    output = np.array([wavelength, inter_flux, new_inter_var]) # put the interpolated data into the new table

    return output # return new table

SN = np.genfromtxt('../../data/spectra/bsnip/sn2002ha-20021102-ui-corrected.flm')

wavelength = SN[:,0]
flux = SN[:,1]
try:
    error = SN[:,2]
except IndexError:
    error = np.array([0])
ivar = genivar(wavelength, flux, error)

newdata = Interpo(wavelength, flux, ivar)

plt.plot(newdata[0,:], newdata[1,:])
plt.show()
