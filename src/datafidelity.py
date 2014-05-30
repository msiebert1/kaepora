############################################################################
#
## This module contains the functions for data fidelity:
## Smoothing functions, clip_data functions, and variance generation
#
## To include these functions, use
#  from datafidelity import *
#
############################################################################

# Import necessary python modules
import numpy as np
from math import *
from scipy import interpolate
import matplotlib.pyplot as plt
import pyfits


############################################################################
#
## Smoothing functions
#
## Function gsmooth() is an inverse variance weighted Gaussian smoothing of spectra
## Optional imputs are smoothing velocity (vexp) and number of sigma (nsig)
## Syntax: new_y_array = gsmooth(x_array, y_array, var_y, vexp = 0.01, nsig = 5.0)

def gsmooth(x_array, y_array, var_y, vexp = 0.005, nsig = 5.0):
    
    # Check for non-zero variance points, and set to 1E-20
    for i in range(len(var_y)):
        if var_y[i] == 0:
            var_y[i] = 1E-20
    
    # Output y-array
    new_y = np.zeros(len(x_array), float)
    
    # Loop over y-array elements
    for i in range(len(x_array)):
        
        # Construct a Gaussian of sigma = vexp*x_array[i]
        gaussian = np.zeros(len(x_array), float)
        sigma = vexp*x_array[i]
        
        # Restrict range to +/- nsig sigma
        sigrange = np.nonzero(abs(x_array-x_array[i]) <= nsig*sigma)
        gaussian[sigrange] = (1/(sigma*sqrt(2*pi)))*np.exp(-0.5*((x_array[sigrange]-x_array[i])/sigma)**2)
        
        # Multiply Gaussian by 1 / variance
        W_lambda = gaussian / var_y
        
        # Perform a weighted sum to give smoothed y value at x_array[i]
        W0 = np.sum(W_lambda)
        W1 = np.sum(W_lambda*y_array)
        new_y[i] = W1/W0

    # Return smoothed y-array
    return new_y

############################################################################
#
# Function clip() tries to identify absorption lines and cosmic rays
# Required input is wavelength, flux and inverse variance arrays

def clip(wave, flux, ivar):
    # Create an array of all ones
    var = np.ones(len(flux), float)
    
    # Create 2 smoothed fluxes, of varying vexp
    sflux = gsmooth(wave, flux, var, 0.002)
    
    # Take the difference of the two fluxes and smooth
    err = abs(flux - sflux)  
    serr = gsmooth(wave, err, var, 0.008)

    # Find the wavelengths that need to be clipped (omitting 5800-6000 region)
    bad_wave = wave[np.where((err/serr > 4) & ((wave < 5800.0) | (wave > 6000.0)))]

    # Find indices for general clipping
    bad = np.array([], int)
    #bad_ranges = [] # don't need it to be a numpy array (A.S.)
    for i in range(len(bad_wave)):
        bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))
        #bad_ranges.append((bad_wave[i]-8, bad_wave[i]+8)) # instead save each point as tuple of desired range (A.S.)

    # Set ivar to 0 for those points and return
    ivar[bad] = 0
    return ivar
    #return bad_ranges # return bad_ranges instead of setting ivar[bad] = 0 (A.S.)

############################################################################
#
# Function to add sky over a wavelength range
#
def addsky(wavelength, flux, error, med_error):

    # Open kecksky spectrum from fits file and create arrays
    sky = pyfits.open('../personal/AdamSnyder/kecksky.fits')
    
    crval = sky[0].header['CRVAL1']
    delta = sky[0].header['CDELT1']
    skyflux = sky[0].data[0]
    start = crval - ceil(0.5*len(skyflux)*delta)
    stop = crval + ceil(0.5*len(skyflux)*delta)
    skywave = [(start+delta*i) for i in range(len(skyflux))]

    # Find wavelength overlap
    good = np.where((wavelength >= skywave[0]) & (wavelength <= skywave[-1]))

    if len(good[0]) == 0:
        return error

    spline_rep = interpolate.splrep(skywave, skyflux)
    add_flux = interpolate.splev(wavelength[good], spline_rep)    

    # Scale sky
    scale = 285*med_error
    add_flux = scale*add_flux

    # Add sky flux to the error
    new_error = error
    new_error[good] = error[good] + add_flux

    return new_error

############################################################################
#
# Function to generate inverse variance for files that are missing this data
#
## Function genivar() generates an inverse variance spectrum.
## Required inputs are an array of wavelengths (wavelength) and an array of corresponding fluxes (flux)
## Optional inputs are velocity of smoothing (vexp) [default 0.0008] and number of sigma (nsig) [default 5.0]
## genivar(wavelength, flux, float vexp = 0.005, float nsig = 5.0)
#

def genivar(wavelength, flux, varflux , vexp = 0.0008, nsig = 5.0): # varflux should be set as the input error?

    # Check to see if it has a variance already
    if len(np.where(varflux == 0)[0]) > 0:
        varflux = np.ones(len(wavelength))
    else:
        ivar = 1 / (varflux**2)
        return ivar
    
    # Smooth original flux
    new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)
    
    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)
    
    # Smooth noise to find the variance
    sm_error = gsmooth(wavelength, error, varflux, vexp, nsig)

    # Test wavelength ranges for kecksky overlap
    test1 = np.where((wavelength >= 5000) & (wavelength <= 6000))
    test2 = np.where((wavelength >= 6000) & (wavelength <= 7000))
    test3 = np.where((wavelength >= 7000) & (wavelength <= 8000))
    if len(test1[0]) > 40:
        med_err = 1.8*np.median(sm_error[test1])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    elif len(test2[0]) > 40:
        med_err = 1.8*np.median(sm_error[test2])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    elif len(test3[0]) > 40:
        med_err = 1.8*np.median(sm_error[test3])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    else:
        sm_error_new = sm_error
        
    # Inverse variance
    ivar = 1/(sm_error_new**2)
    
    # Return generated variance
    return ivar
