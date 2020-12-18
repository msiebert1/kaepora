import numpy as np
from math import *
from scipy import interpolate
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import copy

import scipy.signal as sig

"""This file contains necessary functions for improving raw spectra. Notably, smoothing,
    clipping, variance spectrum generation, and adding extra variance to due sky lines.
"""
def find_vexp(x_array, y_array, var_y=None):
    """ Determines the optimal smoothing parameter (vexp) based on an estimate of SNR
    """
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

    return vexp_auto, SNR

def gsmooth(x_array, y_array, var_y, vexp , nsig = 5.0):
    """Function gsmooth() is an inverse variance weighted Gaussian smoothing of spectra
       Optional imputs are smoothing velocity (vexp) and number of sigma (nsig)
       Syntax: new_y_array = gsmooth(x_array, y_array, var_y, vexp = 0.01, nsig = 5.0)
    """
    
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


def clip(wave, flux, var, vexp, testing=False, filename=None):
    """Uses sigma clipping to find and interpolate over unwanted cosmic rays, and emission lines.
    """
    
    # Create 2 smoothed fluxes, of varying vexp
    sflux = gsmooth(wave, flux, var, vexp)

    # Take the difference of the two fluxes and smooth
    diff = flux - sflux
    err = abs(flux - sflux)
    serr = gsmooth(wave, err, var, 0.008)

    # Avoid clipping absorption at Ca H&K 3934, Na I D 5890/5896, K I 7665/7699, DIB 5780
    regions_avoid = np.where(((wave > 3919.0) & (wave < 3949.0)) | ((wave > 5875.0) & (wave < 5911.0)) | 
                             ((wave > 7650.0) & (wave < 7714.0)) | ((wave > 5765.0) & (wave < 5795.0)))
    tol1 = 5.
    tol2 = 10.
    bad_inds_pos_avoid = np.where((diff[regions_avoid] > 0) & (err[regions_avoid]/serr[regions_avoid] > tol1))
    bad_inds_neg_avoid = np.where((diff[regions_avoid] < 0) & (err[regions_avoid]/serr[regions_avoid] > tol2))

    regions_normal = np.where((wave <= 3919.0) | ((wave >= 3949.0) & (wave <= 5765.0)) | 
                              ((wave >= 5795.0) & (wave <= 5875.0)) | ((wave >= 5911.0) & (wave <= 7650.0)) | 
                              (wave >= 7714.0))

    bad_inds_normal = np.where((err[regions_normal]/serr[regions_normal] > tol1))

    bad_wave_normal = wave[regions_normal][bad_inds_normal]
    bad_wave_pos_avoid = wave[regions_avoid][bad_inds_pos_avoid]
    bad_wave_neg_avoid = wave[regions_avoid][bad_inds_neg_avoid]

    # Find indices for general clipping
    bad_ranges = [] 
    buff = 8
    for i in range(len(bad_wave_normal)):
        if bad_wave_normal[i] + buff < wave[-1] and bad_wave_normal[i] - buff >= wave[0]:
            bad_ranges.append((bad_wave_normal[i]-buff, bad_wave_normal[i]+buff))
        # elif bad_wave_normal[i] + buff >= wave[-1]:
        #     bad_ranges.append((bad_wave_normal[i]-buff, wave[-2]))
        # elif bad_wave_normal[i] - buff < wave[0]:
        #     bad_ranges.append((wave[0], bad_wave_normal[i]+buff))


    for i in range(len(bad_wave_pos_avoid)):
        bad_ranges.append((bad_wave_pos_avoid[i]-buff, bad_wave_pos_avoid[i]+buff))
    for i in range(len(bad_wave_neg_avoid)):
        bad_ranges.append((bad_wave_neg_avoid[i]-buff, bad_wave_neg_avoid[i]+buff))

    # print bad_ranges
    if testing:
        scale = 10./np.amax(flux)
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
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

    for i, wave_tuple in enumerate(bad_ranges):
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
        # plt.savefig('../../Foundation/Spectra/Clipping_Plots/'+filename+'_clipping.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()

    return wave, flux, var # return bad_ranges instead of setting ivar[bad] = 0 (A.S.)
    

def clipmore(wave, flux, ivar) :
    """ This function is no longer used, Na and other specific regions are 
        handled by the clip() function
    """

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
        

def addsky(wavelength, flux, error, med_error, source = None):

    # Open kecksky spectrum from fits file and create arrays
    # sky = pyfits.open('../personal/AdamSnyder/kecksky.fits')
    if source.lower() == 'bsnip':
        sky = pyfits.open('../data/sky_spectra/licksky.fits')
    else:
        sky = pyfits.open('../data/sky_spectra/kecksky.fits')

    crval = sky[0].header['CRVAL1']
    delta = sky[0].header['CDELT1']
    if source.lower() == 'bsnip':
        skyflux = sky[0].data
    else:
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

    scale = med_error/np.median(add_flux)

    add_flux = 0.004*scale*add_flux

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

def genivar(wavelength, flux, varflux, vexp = 0.002, nsig = 5.0, testing=False, sky=True, source=None):
    """ The primary function for variance spectrum generation. Spectrum is smoothed, residual spectrum 
        is generated, residual spectrum is smoothed, sky noise is added in, error spectrum is scaled
        by a linear function of wavelength determined from the CfA data, and the ivar spectrum and SNR 
        are returned.
    """

    if varflux is not None and not testing:
        error = np.sqrt(varflux)
        SNR = np.median(flux / error)
        ivar = 1. / (varflux)
        return ivar, SNR

    # Smooth original flux

    new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)
    
    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)
    
    # Smooth noise to find the variance
    sm_error = gsmooth(wavelength, error, varflux, .015, nsig)

    # Test wavelength ranges for kecksky overlap
    test1 = np.where((wavelength >= 5000) & (wavelength <= 6000))
    test2 = np.where((wavelength >= 6000) & (wavelength <= 7000))
    test3 = np.where((wavelength >= 7000) & (wavelength <= 8000))

    if sky:
        if len(test1[0]) > 40:
            med_err = np.median(sm_error[test1])
            sm_error_new = addsky(wavelength, flux, sm_error, med_err, source=source)
        elif len(test2[0]) > 40:
            med_err = np.median(sm_error[test2])
            sm_error_new = addsky(wavelength, flux, sm_error, med_err, source=source)
        elif len(test3[0]) > 40:
            med_err = np.median(sm_error[test3])
            sm_error_new = addsky(wavelength, flux, sm_error, med_err, source=source)
        else:
            sm_error_new = sm_error
    else:
        sm_error_new = sm_error

    if varflux is not None:
        real_error = np.sqrt(varflux)
        sm_real_error = gsmooth(wavelength, real_error, None, .008)
        error_scales = sm_real_error/sm_error
        scale_fit = np.polyfit(wavelength, error_scales, 1)
        scale_func = scale_fit[0]*(wavelength) + scale_fit[1]
    

    scale_func = -5.01307400136e-05*(wavelength) + 1.61999679746 #RESULTING COEFFS FROM RECENT FITS
    scaled_error = scale_func*sm_error_new

    if testing:
        scale = 1.
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
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
        plt.plot(wavelength, scale*sm_error_new, linewidth = 2, color = 'orange', label='Smoothed Absolute Residuals')
        plt.plot(wavelength, scale*scale_func*sm_error_new, linewidth = 2, color = 'magenta', label='Smoothed Absolute Residuals Corrected')
        if varflux is not None:
            plt.plot(wavelength, scale*real_error, linewidth = 2, color = 'g', label='Real Uncertainty')
        plt.ylabel('Uncertainty (Scaled Flux)', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([wavelength[0]-200,wavelength[-1]+200])
        plt.legend(loc=1, fontsize=20)
        # plt.savefig('../../../Paper_Drafts/reprocessing_updated/genvar.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()
        if varflux is not None:
            plt.plot(wavelength, error_scales)
            plt.plot(wavelength, scale_func)
            plt.show()
            
    sm_var_new = scaled_error*scaled_error

    # Inverse variance
    ivar = 1./(sm_var_new)
    SNR = np.median(new_flux / scaled_error)
    return ivar, SNR
