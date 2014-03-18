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
from math import floor,ceil
#import sqlite3 as sq3
#import msgpack

"""
README:

This part of code is written for composite spectra preparation.
It does deredding, deredshifting, and interpolation.
Using the function prep:
INPUT :
SPECTRUM : table containing the original wavelength and flux
FILE_NAME : containing name of supernova


OUTPUT:
NEW_DATA: table containing the processed wavelength, flux and variance

"""

def ReadParam():
    #Read in : table containing sn names, redshifts, etc.
    sn_parameter = np.genfromtxt('../data/cfa/cfasnIa_param.dat',dtype = None)

    return sn_parameter

def ReadExtin():
    #table containing B and V values for determining extinction -> dereddening due to milky way
    sne = np.genfromtxt('extinction.dat', dtype = None)

    return sne



"""
NOTE:
Using IRSA
"""

"""
Note: using parameters.dat file which was created from paramaters.py
parameters.py is designed to pull all relevant parameters for SNe spectra from online databases
via astroquery and place it all in a table to be pulled from later;
it would ideally do that following:
-read in all files we want to prep
-truncate file name to format "SNyear"(this is how astroquery searches for SN)
-get relevent data into table, following a format like:
SN name		Host Galaxy		Redshift	B	V	De-redden factor	CarbonPos/Neg

####
NOTE:Currently only has SN_name, B, and V values for purposes of Dereddening due to Milky way dust
####
"""

#deredden spectra to milky way
#deredshift the spectra
#deredden to host galaxy

def dered(sn_param,sne,filename,wave,flux):
    for j in range(len(sn_param)):#go through list of SN parameters
        sn = sn_param[j][0] #get relevent parameters needed for calculations
        z = sn_param[j][1]  # redshift value
        if sn in filename:#SN with parameter matches the path
            b = sne[j][1].astype(float)
            v = sne[j][2].astype(float)
            bv = b-v
            r = v/bv
#            print "B(%s)-V(%s)=%s"%(b,v,bv)
#            print "R(v) =",r
            #or use fm07 model
            #test1 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='ccm89')
            #test2 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='od94')
            flux *= ex.reddening(wave,ebv = bv, r_v = 3.1, model='f99')
            wave /= (1+z)

            #print "de-reddened by host galaxy\n",flux*ex.reddening(wave,ebv = 0, r_v = r, model='f99')
            #host *= ex.reddening(wave,ebv = bv, r_v = r, model='f99')

    return [wave,flux]


##############################################################################################################################################
##############################################################################################################################################

# Data Interpolation

"""
NOTE:
This function inputs three lists containing wavelength, flux and variance.
The output will be a Table with all the fitted values.
You can change the mininum and maximum of wavelength to output, as well as pixel size in the first few lines.
For the spectra that does not cover the whole range of specified wavelength,
we output the outside values as NAN
"""


def Interpo (wave, flux, variance) :
    wave_min = 1500
    wave_max = 12000
    dw = 2

    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)
    wavelength = np.arange(ceil(wave_min), floor(wave_max), dtype=int, step=dw) #creates N equally spaced wavelength values
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

    missing_data = np.where((wavelength < lower) | (wavelength > upper))
    inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
    inter_var[missing_data] =  float('NaN')

    output = np.array([wavelength, inter_flux, inter_var]) # put the interpolated data into the new table

    return output # return new table


    # Get the Noise for each spectra

def getsnr(flux, variance) :

    snr = flux/variance
    snr_med = np.median(snr)
    return snr_med

from datafidelity import *  # Get variance from the datafidelity outcome



def compprep(spectrum,file_name):
    sn_parameter = ReadParam()
    sne = ReadExtin()
    newdata = []

    old_wave = spectrum[:,0]	    #wavelengths
    old_flux = spectrum[:,1] 	#fluxes
    #old_var  = spectrum[:,2]  #errors

    old_var = genvar(old_wave, old_flux) #variance
    snr = getsnr(old_flux, old_var)
    print 'S/N ratio', file_name, snr

    new_spectrum = dered(sn_parameter, sne, file_name, old_wave, old_flux)
    new_wave = new_spectrum[0]
    new_flux = new_spectrum[1]
    new_var  = genvar(new_wave, new_flux) #variance
    #var = new_flux*0+1
    newdata = Interpo(new_wave, new_flux, new_var) # Do the interpolation
#    print 'new spectra',newdata
    return newdata, snr

