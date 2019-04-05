import os
import glob
from specutils import extinction as ex
from specutils import Spectrum1D
from astropy import units as u
import test_dered
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import math
import scipy.optimize as opt
import copy

""" This file contains various functions for homogenizing our dataset. 
    The compprep() function is called before a spectrum is added to the
    spectral table in kaepora. This takes the raw spectra and homogenizes 
    them via the methods outlined in Siebert et al. 2019.
"""

def ReadParam():
    #Read in : table containing sn names, redshifts, etc.
    sn_param = np.genfromtxt('../data/cfa/cfasnIa_param.dat', dtype=None)
    sn = []
    z = []
    for i in range(len(sn_param)):
        sn.append(sn_param[i][0])  # get relevent parameters needed for calculations
        z.append(sn_param[i][1])  # redshift value
    return z


def ReadExtin(file):
    #table containing B and V values for determining extinction -> dereddening due to milky way
    sne = np.genfromtxt(file, dtype=None)

    return sne


def dered(sne, snname, wave, flux):
    """This function is deprecated. compprep() now uses 
        test_dered.dered() (yes that is a dumb file name)
    """
    for j in range(len(sne)):  # go through list of SN parameters
        sn = sne[j][0]
        if sn in snname:  # SN with parameter matches the path
            b = sne[j][1].astype(float)
            v = sne[j][2].astype(float)
            bv = b-v
#            print "B(%s)-V(%s)=%s"%(b,v,bv)
#            print "R(v) =",r
            #or use fm07 model
            #test1 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='ccm89')
            #test2 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='od94')
            flux *= ex.reddening(wave, ebv=bv, r_v=3.1, model='f99')
#            wave /= (1+z)

            #print "de-reddened by host galaxy\n",flux*ex.reddening(wave,ebv = 0, r_v = r, model='f99')
            #host *= ex.reddening(wave,ebv = bv, r_v = r, model='f99')

    return flux

def host_correction(sne, snname, wave, flux):
    """This function is deprecated. composite.py now uses 
        test_dered.host_correction() (yes that is a dumb file name)
    """
    for j in range(len(sne)):  # go through list of SN parameters
        sn = sne[j][0]
        if sn in snname:  # SN with parameter matches the path
            a_v = sne[j][2].astype(float)
            flux *= ex.reddening(wave, ebv=3.1*a_v, r_v=3.1, model='f99')
    return flux


# Data Interpolation

import datafidelity as df

def Interpo_flux_conserving(wave, flux, ivar, dw=2, testing=False):
    """This is a an interpolation algorithm that does trapezoidal integration
        to conserve flux. The variance is then propagated correctly. Since 
        interpolation always introduces correlation in the variance spectrum,
        we ignore  this correlation byscale the original variance spectrum 
        to the new variance spectrum.
    """
    var = 1./ivar
    pixel_scale = np.median(np.diff(wave))
    # print pixel_scale

    wave_min = 999
    wave_max = 12001
    wavelength_mids = np.arange(math.ceil(wave_min), math.floor(wave_max),
                           dtype=int, step=dw)
    lower = wave[0]
    upper = wave[-1]

    good_data = np.where((wave >= lower) & (wave <= upper))
    influx = inter.splrep(wave[good_data], flux[good_data])
    invar = inter.splrep(wave[good_data], var[good_data])

    low_int = math.ceil(lower)
    up_int = math.ceil(upper)
    wave_final = []
    flux_final = []
    var_final = []

    for mid_point in wavelength_mids:
        inter_flux_mid_left = float(inter.splev(mid_point, influx, ext = 3))
        inter_var_mid_left = float(inter.splev(mid_point, invar, ext = 3))
        inter_flux_mid_right = float(inter.splev(mid_point+dw, influx, ext = 3))
        inter_var_mid_right = float(inter.splev(mid_point+dw, invar, ext = 3))
        waves_to_sum = np.where((wave > mid_point) & (wave < mid_point+dw))[0]
        wave_final.append((mid_point + mid_point +dw)/2.)
        if (mid_point >= lower and (mid_point+dw) <= upper):
            waves = [mid_point]
            fluxes = [inter_flux_mid_left]
            variances = [inter_var_mid_left]
            for i, w in enumerate(waves_to_sum):
                waves.append(wave[w])
                fluxes.append(flux[w])
                variances.append(var[w])
            waves.append(mid_point+dw)
            fluxes.append(inter_flux_mid_right)
            variances.append(inter_var_mid_right)
            new_point = np.trapz(fluxes, x=waves)
            flux_final.append(new_point)
            diffs = np.diff(waves)
            var_tot = 0.
            for i in range(len(variances)-1):
                v1 = variances[i]
                v2 = variances[i+1]
                var_tot = var_tot + (diffs[i]**2.)*(v1+v2)
            var_tot = var_tot*.25
            var_final.append(var_tot)
        else:
            flux_final.append(0.)
            var_final.append(0.)

    inter_var = inter.splev(wave_final, invar, ext = 3)
    s = scale_composites_in_range(var_final, inter_var)[0]
    scale = 1./s
    var_uncorrelated = scale*inter_var
    ivar_final = 1./var_uncorrelated

    wave_final = np.asarray(wave_final)
    flux_final = np.asarray(flux_final)
    ivar_final = np.asarray(ivar_final)
    if testing:
        var_final = np.asarray(var_final)

    missing_data = np.where((wave_final < lower+2) | (wave_final > upper-2))# padding to avoid interpolation errors near edges
    flux_final[missing_data] = float('NaN')
    ivar_final[missing_data] = float('NaN')
    if testing:
        var_final[missing_data] = float('NaN')

    output = np.array([wave_final, flux_final, ivar_final])

    if testing:
        interp_wave = output[0,:]
        interp_flux = output[1,:]
        interp_ivar = output[2,:]
        # print scale
        plt.plot(wave,var)
        plt.plot(interp_wave,var_final)
        plt.plot(interp_wave,1./interp_ivar)
        plt.xlim([7000,7100])
        plt.ylim([-.05e-32,2.e-32])
        plt.show()

    return output, scale, var_final


def Interpo (wave, flux, ivar):
    """ This is no longer used. Does a simple interpolation of the flux and variance.
    This does not conserve flux or propagate variance correctly.
    """
    wave_min = 1000
    wave_max = 12000
    dw = 2

    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)
    wavelength = np.arange(math.ceil(wave_min), math.floor(wave_max),
                           dtype=int, step=dw)  # creates N equally spaced wavelength values
    inter_flux = []
    inter_ivar = []
    output = []

    lower = wave[0]  # Find the area where interpolation is valid
    upper = wave[-1]

    good_data = np.where((wave >= lower) & (wave <= upper))  #creates an array of wavelength values between minimum and maximum wavelengths from new spectrum

    influx = inter.splrep(wave[good_data], flux[good_data])  # creates b-spline from new spectrum

    inivar = inter.splrep(wave[good_data], ivar[good_data])  # doing the same with the inverse varinces

    # extrapolating returns edge values
    inter_flux = inter.splev(wavelength, influx, ext = 3)	 # fits b-spline over wavelength range
    inter_ivar = inter.splev(wavelength, inivar, ext = 3)   # doing the same with errors

    inter_ivar[inter_ivar < 0] = 0  

    missing_data = np.where((wavelength < lower) | (wavelength > upper))
    inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
    inter_ivar[missing_data] = float('NaN')

    output = np.array([wavelength, inter_flux, inter_ivar])  # put the interpolated data into the new table

    return output  # return new table



def getsnr(flux, ivar):
    """Returns the approximate SNR of a spectrum given its flux and ivar data
    """
    sqvar = map(math.sqrt, ivar)
    snr = flux/(np.divide(1.0, sqvar))
    snr_med = np.median(snr)
    return snr_med

def scale_composites_in_range(data, comp):
    """Finds the scale factor the minimizes the difference between two 
        spectra (data and comp)
    """
    scales = []
    guess = 1.
    s = opt.minimize(sq_residuals_in_range, guess, args = (data, comp), 
                 method = 'Nelder-Mead').x
    return s

def sq_residuals_in_range(s, data, comp):
    """Calculates the sum of the square residuals between arrays data and comp
    """
    data = s*data
    res = data - comp
    sq_res = res*res
    return np.sum(sq_res)


def compprep(spectrum, sn_name, z, source, use_old_error=True, testing=False):
    """ Performs clipping, deredshifting, variance spectrum generation, MW extinction correction,
        and interpolation. If testing is True, several plots will be made to assess the quality 
        of this processing.
    """
    old_wave = spectrum[:, 0]	    # wavelengths
    old_flux = spectrum[:, 1] 	# fluxes
    try:
        old_error = spectrum[:, 2]  # check if supernovae has error array
    except IndexError:
        old_error = None  # if not, set default
    if sn_name == '2011fe' and source == 'other':
        old_error = np.sqrt(old_error)
    if old_error is not None:
        old_var = old_error**2.
    else:
        old_var = None

    if old_var is not None:
        num_non_zeros = np.count_nonzero(old_var)
        if len(old_var) - num_non_zeros > 100:
            old_var = None
        elif old_var[-1] == 0.:
            old_var[-1] = old_var[-2]
        elif True in np.isnan(old_var):
            nan_inds = np.transpose(np.argwhere(np.isnan(old_var)))[0]
            for ind in nan_inds:
                if ind != 0:
                    old_var[ind] = old_var[ind-1]
                else:
                    old_var[ind] = old_var[ind+1]

    # if testing:
    #     plt.plot(old_wave, old_flux)
    #     plt.plot(old_wave/(1.+z), old_flux)
    #     plt.plot(old_wave*(1.+z), old_flux)
    #     plt.xlim(5800,6000)
    #     # plt.show()
    #     if old_var is not None:
    #         plt.plot(old_wave, old_var)
    #         plt.show()
    # old_var = None
    vexp, SNR = df.find_vexp(old_wave, old_flux, var_y=old_var)
    if testing:
        print vexp, SNR

    if source != 'csp': #already deredshifted
        old_wave = old_wave/(1.+z) #deredshift for clipping 
        
    old_wave, old_flux, old_var = df.clip(old_wave, old_flux, old_var, vexp, testing=testing) #clip emission/absorption lines
    old_wave = old_wave*(1.+z) #reredshift for MW extinction correction 
    temp_ivar, SNR = df.genivar(old_wave, old_flux, old_var, vexp=vexp, testing=testing, source=source)  # generate inverse variance

    if testing:
        print SNR

    if old_var is not None:
        old_ivar = 1./old_var
    else:
        old_ivar = temp_ivar
    # snr = getsnr(old_flux, old_ivar)

    if source == 'cfa':  # choosing source dataset
#        z = ReadParam()
        sne = ReadExtin('extinction.dat')
    if source == 'bsnip':
        sne = ReadExtin('extinctionbsnip.dat')
    if source == 'csp':
        sne = ReadExtin('extinctioncsp.dat')
    if source == 'uv':
        sne = ReadExtin('extinctionuv.dat')
    if source == 'other':
        sne = ReadExtin('extinctionother.dat')
    if source == 'swift_uv':
        sne = ReadExtin('extinctionswiftuv.dat')
    if source == 'foley_hst':
        sne = ReadExtin('extinctionhst.dat')

#     host_reddened = ReadExtin('../data/info_files/ryan_av.txt')
    newdata = []
    old_wave = old_wave*u.Angstrom        # wavelengths
    old_flux = old_flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(old_wave, old_flux)
    spec1d_ivar = Spectrum1D.from_array(old_wave, old_ivar)
    dered_flux, dered_ivar = test_dered.dered(sne, sn_name, spec1d.wavelength, spec1d.flux, spec1d_ivar.flux)  # Dereddening (see if sne in extinction files match the SN name)
#     new_flux = host_correction(sne, sn_name, old_wave, new_flux)

    # new_flux = old_flux

    if testing:
        new_flux_plot = copy.deepcopy(dered_flux)
        new_ivar_plot = copy.deepcopy(dered_ivar)
        old_wave_plot = copy.deepcopy(old_wave)

    new_flux = dered_flux.value
    new_ivar = dered_ivar.value
    old_wave = old_wave.value

    if testing:
        av_specific = 0.2384 #2005lz
        av_specific = 0.4089 #2007af
        r_v=2.5
        new_flux_host, new_ivar_host = test_dered.host_correction(av_specific, r_v, sn_name, old_wave_plot, new_flux_plot, new_ivar_plot)
        new_flux_host = new_flux_host.value
        old_flux = old_flux.value

        s = scale_composites_in_range(new_flux, old_flux)
        new_flux_scaled = s*new_flux
        s = scale_composites_in_range(new_flux_host, old_flux)
        new_flux_host_scaled = s*new_flux_host

        valid_data = np.where(old_wave > 4000)
        norm = 10./np.nanmax(new_flux_host_scaled[valid_data])
        old_flux_norm = old_flux*norm
        new_flux_norm = new_flux_scaled*norm
        new_flux_host_norm = new_flux_host_scaled*norm

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
        plt.plot(old_wave, old_flux_norm, linewidth = 2, color = '#000080', label='Before Dereddening')
        plt.plot(old_wave, new_flux_norm, linewidth = 2, color = 'gold', label='Milky Way Corrected')
        plt.plot(old_wave, new_flux_host_norm, linewidth = 2, color = '#d95f02', label='Host Corrected')
        plt.ylabel('Relative Flux', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([old_wave[0]-200,old_wave[-1]+200])
        plt.legend(loc=1, fontsize=20)
        # plt.savefig('../../../Paper_Drafts/reprocessing_updated/red_corr.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()
        # plt.plot(old_wave, old_ivar)
        # plt.plot(old_wave, new_ivar)
        # plt.show()

    new_wave = old_wave/(1.+z)  # Deredshifting

    if not use_old_error:
        new_var = None
    else:
        new_var = old_var  # Placeholder if it needs to be changed
    #var = new_flux*0+1
    # newdata = Interpo(new_wave, new_flux, new_ivar)  # Do the interpolation
    newdata, scale, var_final = Interpo_flux_conserving(new_wave, new_flux, new_ivar, testing=testing)

    if testing:
        # newdata_test = Interpo(new_wave, new_flux_host_norm, new_ivar)
        newdata_test, scale, var_final = Interpo_flux_conserving(new_wave, new_flux_host_norm, new_ivar)
        interp_wave = newdata_test[0,:]
        interp_flux = newdata_test[1,:]
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
        plt.plot(new_wave, 2.*new_flux_host_norm, linewidth = 2, color = '#d95f02', label='Before Interpolation')
        plt.plot(interp_wave, interp_flux, linewidth = 2, color = 'darkgreen', label='After Interpolation')
        plt.ylabel('Relative Flux', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([new_wave[0]-200,new_wave[-1]+200])
        plt.legend(loc=1, fontsize=20)
        # plt.savefig('../../../Paper_Drafts/reprocessing_updated/interp_deredshift.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()

    return newdata, SNR
