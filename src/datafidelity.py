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
# import pyfits
import astropy.io.fits as pyfits
import copy

import scipy.signal as sig


############################################################################
#
## Smoothing functions
#
## Function gsmooth() is an inverse variance weighted Gaussian smoothing of spectra
## Optional imputs are smoothing velocity (vexp) and number of sigma (nsig)
## Syntax: new_y_array = gsmooth(x_array, y_array, var_y, vexp = 0.01, nsig = 5.0)

def find_vexp(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = gsmooth(x_array, y_array, var_y, .002)
        SNR = np.median(new_y_init / error)
    else:
        new_y_init = gsmooth(x_array, y_array, var_y, .002) #this smoothing should get in right ballpark
        error = np.absolute(y_array - new_y_init)
        sm_error = gsmooth(x_array, error, var_y, .008)
        SNR = np.median(new_y_init / sm_error)

    #TODO: interpolate a function of SNR
    # vexp_line = np.polyfit([2.5, 80], [.0045, .001], 1)
    # coeff_0 = vexp_line[0]
    # coeff_1 = vexp_line[1]
    # results from above:
    coeff_0 = -4.51612903e-05
    coeff_1 = 4.61290323e-03
    vexp_auto = coeff_0*SNR + coeff_1

    if SNR < 2.5:
        vexp_auto = .0045
    if SNR > 80:
        vexp_auto = .001
    # if SNR < 5:
    #     vexp_auto = .0045 
    # elif 5 <= SNR < 20:
    #     vexp_auto = .004
    # elif 20 <= SNR < 40:
    #     vexp_auto = .003
    # elif 40 <= SNR < 60:
    #     vexp_auto = .002
    # elif 60 <= SNR < 100:
    #     vexp_auto = .0015
    # else:
    #     vexp_auto = .001

    return vexp_auto, SNR

def gsmooth(x_array, y_array, var_y, vexp , nsig = 5.0):
    
    # Check for zero variance points, and set to 1E-20

    if var_y is None:
        var_y = 1.e-31*np.ones(len(y_array))
    
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

def clip(wave, flux, var, vexp, testing=False):
    # Create an array of all ones
    # var = np.zeros(len(flux), float)
    
    # Create 2 smoothed fluxes, of varying vexp
    sflux = gsmooth(wave, flux, var, vexp)
    # print 'sflux', sflux
    # Take the difference of the two fluxes and smooth
    err = abs(flux - sflux)
    serr = gsmooth(wave, err, var, 0.008)

    # Find the wavelengths that need to be clipped (omitting 5800-6000 region)
    tol1 = 5.5
    tol2 = 5.5
    tol_na1 = 10.5
    index = np.where(((err/serr > tol1) & (wave < 5800.0)) | ((err/serr > tol2) & (wave > 6000.0)))
    index_na = np.where((err/serr > tol_na1) & ((wave > 5800.0) & (wave < 6000.0)))
    # index = np.where(((flux/sflux > tol1*serr) & (wave < 5800.0)) | ((flux/sflux > tol2*serr) & (wave > 6000.0)))
    # index_na = np.where((flux/sflux > tol_na1*serr) & ((wave > 5800.0) & (wave < 6000.0)))
    bad_wave = wave[index]
    bad_wave_na = wave[index_na]

    # Find indices for general clipping
    # bad = np.array([], int)
    bad_ranges = [] # if don't need it to be a numpy array (A.S.)
    buff = 8
    for i in range(len(bad_wave)):
        # bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))
        if bad_wave[i] + buff < wave[-1] and bad_wave[i] - buff >= wave[0]:
            bad_ranges.append((bad_wave[i]-buff, bad_wave[i]+buff))
        elif bad_wave[i] + buff >= wave[-1]:
            bad_ranges.append((bad_wave[i]-buff, wave[-2]))
        elif bad_wave[i] - buff < wave[0]:
            bad_ranges.append((wave[0], bad_wave[i]+buff))



    for i in range(len(bad_wave_na)):
        # bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))
        bad_ranges.append((bad_wave_na[i]-buff, bad_wave_na[i]+buff))

    # Set ivar to 0 for those points and return
#    ivar[bad] = 0
#    return ivar
    if testing:
        scale = 10./np.amax(flux)
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
        # ax.xaxis.set_ticks(np.arange(np.round(wave[0],-3),np.round(wave[-1],-3),1000))
        plt.yticks(fontsize = 20)
        plt.tick_params(
            which='major', 
            bottom='on', 
            top='on',
            left='on',
            right='on',
            length=10)
        plt.tick_params(
            which='minor', 
            bottom='on', 
            top='on',
            left='on',
            right='on',
            length=5)
        plt.plot(wave, scale*flux, 'r-', linewidth = 2, label='Before Clipping')

    for wave_tuple in bad_ranges:
#        print wave_tuple
        clip_points = np.where((wave > wave_tuple[0]) & (wave < wave_tuple[1]))
        #make sure not at edge of spectrum
        flux[clip_points] = np.interp(wave[clip_points], [wave_tuple[0], wave_tuple[1]], [flux[clip_points[0][0]-1], flux[clip_points[0][-1]+1]])

        #deweight data (but not to 0), somewhat arbitrary
        if var is not None:
            var[clip_points] = np.interp(wave[clip_points], [wave_tuple[0], wave_tuple[1]], [var[clip_points[0][0]-1], var[clip_points[0][-1]+1]])

    if testing:
        plt.plot(wave, scale*flux, linewidth = 2, color = '#000080', label='After Clipping')
        plt.plot(wave, scale*sflux, color='gold', linewidth = 2, label='After Smoothing')
        plt.ylabel('Relative Flux', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([wave[0]-200,wave[-1]+200])
        plt.ylim([0,1.05*10])
        plt.legend(loc=1, fontsize=20)
        # plt.savefig('../../../Paper_Drafts/reprocessing_updated/clipping.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()

    return wave, flux, var # return bad_ranges instead of setting ivar[bad] = 0 (A.S.)
    
# clip function specifically for lines that need ivar to set to nearby ones
# Currently only for Na lines (5800-5900)
def clipmore(wave, flux, ivar) :

    ind = np.where((wave > 5800.0 ) & (wave < 6000.0 )) # what region to look at 
#    print ind[0]
    if len(ind[0]) == 0 : # Doesn't have spectra at this range
        print 'no data'
        return ivar
    else: 
        wmin = ind[0][0]
        wmax = ind[0][-1]
#        print wmin, wmax
        # Create an array of all ones
        var = np.ones(len(flux), float)
    
        # Create 2 smoothed fluxes, of varying vexp
        sflux = gsmooth(wave, flux, var, 0.002)
    
        # Take the difference of the two fluxes and smooth
        err = abs(flux - sflux)  
#        serr = gsmooth(wave, err, var, 0.008)
        
        mederr = np.median(err[ind])
        if mederr/ err[wmin-1] >2 or mederr / err[wmax+1] > 2 :
#        medivar = np.median(ivar[ind])
#        if ivar[wmin-1] / medivar > 4 or ivar[wmax+1] / medivar > 4 : # criterion: mederr/outerr > 2
            print 'sodium clipped!'
#            print ivar[ind]
            ivar[ind] = mederr**-2.0
#            ivar[ind] = medivar
#            print ivar[ind]

        return ivar
        
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
    start = crval
    stop = crval + ceil(len(skyflux)*delta)
    skywave = [(start+delta*i) for i in range(len(skyflux))]
    # Find wavelength overlap
    good = np.where((wavelength >= skywave[0]) & (wavelength <= skywave[-1]))

    if len(good[0]) == 0:
        return error

    spline_rep = interpolate.splrep(skywave, skyflux)
    add_flux = interpolate.splev(wavelength[good], spline_rep)    

    # Scale sky
    # scale = 285*med_error
    # print med_error
    scale = 50.*med_error #fudge factor provides reasonable scaling of sky lines
    add_flux = scale*add_flux

    # Add sky flux to the error
    new_error = copy.deepcopy(error)
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

def genivar(wavelength, flux, varflux, vexp = 0.002, nsig = 5.0, testing=False):
    # Check to see if it has a variance already

    #UNCOMMENT WHEN DONE TESTING
    if varflux is not None and not testing:
        error = np.sqrt(varflux)
        SNR = np.median(flux / error)
        ivar = 1. / (varflux)
        return ivar, SNR

    # Smooth original flux
    # new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)

    new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)
    
    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)
    
    # Smooth noise to find the variance
    # sm_error = gsmooth(wavelength, error, varflux, vexp, nsig)
    sm_error = gsmooth(wavelength, error, varflux, .015, nsig)

    # Test wavelength ranges for kecksky overlap
    test1 = np.where((wavelength >= 5000) & (wavelength <= 6000))
    test2 = np.where((wavelength >= 6000) & (wavelength <= 7000))
    test3 = np.where((wavelength >= 7000) & (wavelength <= 8000))
    if len(test1[0]) > 40:
        med_err = np.median(sm_error[test1])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    elif len(test2[0]) > 40:
        med_err = np.median(sm_error[test2])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    elif len(test3[0]) > 40:
        med_err = np.median(sm_error[test3])
        sm_error_new = addsky(wavelength, flux, sm_error, med_err)
    else:
        sm_error_new = sm_error

    if varflux is not None:
        real_error = np.sqrt(varflux)
        sm_real_error = gsmooth(wavelength, real_error, None, .008)
        error_scales = sm_real_error/sm_error
        scale_fit = np.polyfit(wavelength, error_scales, 1)
        # print scale_fit
        scale_func = scale_fit[0]*(wavelength) + scale_fit[1]
    

    scale_func = -5.01307400136e-05*(wavelength) + 1.61999679746 #RESULTING COEFFS FROM RECENT FITS
    scaled_error = scale_func*sm_error_new

    if testing:
        scale = 10./np.amax(flux)
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
        # ax.xaxis.set_ticks(np.arange(np.round(wave[0],-3),np.round(wave[-1],-3),1000))
        plt.yticks(fontsize = 20)
        plt.tick_params(
            which='major', 
            bottom='on', 
            top='on',
            left='on',
            right='on',
            length=10)
        plt.tick_params(
            which='minor', 
            bottom='on', 
            top='on',
            left='on',
            right='on',
            length=5)

        plt.plot(wavelength, scale*error, linewidth = 2, color = '#000080', label='Absolute Residuals')
        # if varflux is not None:
        #     plt.plot(wavelength, sm_real_error, linewidth = 2, color = 'g', label='Real Error')
        plt.plot(wavelength, scale*sm_error_new, linewidth = 2, color = 'orange', label='Smoothed Absolute Residuals')
        plt.plot(wavelength, scale*scale_func*sm_error_new, linewidth = 2, color = 'magenta', label='Smoothed Absolute Residuals Corrected')
        if varflux is not None:
            plt.plot(wavelength, scale*real_error, linewidth = 2, color = 'g', label='Real Uncertainty')
        plt.ylabel('Uncertainty (Scaled Flux)', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([wavelength[0]-200,wavelength[-1]+200])
        plt.legend(loc=1, fontsize=20)
        plt.savefig('../../../Paper_Drafts/reprocessing_updated/genvar.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()
        if varflux is not None:
            plt.plot(wavelength, error_scales)
            plt.plot(wavelength, scale_func)
            plt.show()
            
    sm_var_new = scaled_error*scaled_error

    # Inverse variance
    ivar = 1./(sm_var_new)
    SNR = np.median(new_flux / scaled_error)
    # Return generated variance
    return ivar, SNR
