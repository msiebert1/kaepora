#!/usr/bin/env python
import datafidelity as df
import matplotlib.pyplot as plt
import numpy as np 
import copy
from scipy.integrate import simps
import random 
# import pyphot
import kaepora as kpora
import copy
from specutils import extinction as ex
from specutils import Spectrum1D
from astropy import units as u
from optparse import OptionParser
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import json
import math
import matplotlib

def autosmooth(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = df.gsmooth(x_array, y_array, var_y, .001)
        SNR = np.mean(y_array / error)
        # SNR = np.mean(new_y_init / error)
    else:
        var_y = np.ones(len(x_array))
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002)
        error = np.absolute(y_array - new_y_init)
        sm_error = df.gsmooth(x_array, error, var_y, .008)
        SNR = np.median(new_y_init / sm_error)

    if SNR < 5:
        vexp_auto = .0045 #temp value, need to fine tune using SNR in excel spreadsheet
    elif 5 <= SNR < 20:
        vexp_auto = .004
    elif 20 <= SNR < 40:
        vexp_auto = .003
    elif 40 <= SNR < 60:
        vexp_auto = .002
    elif 60 <= SNR < 100:
        vexp_auto = .0015
    else:
        vexp_auto = .001

    return vexp_auto, SNR

def find_vexp(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002)
        SNR = np.median(new_y_init / error)
    else:
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002) #this smoothing should get in right ballpark
        error = np.absolute(y_array - new_y_init)
        sm_error = df.gsmooth(x_array, error, var_y, .008)
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

def gsmooth(x_array, y_array, var_y=None, vexp=.002):
    sm_flux = df.gsmooth(x_array, y_array, var_y, vexp)
    return sm_flux

def find_vexp_ryan(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002)
        SNR = np.median(new_y_init / error)
    else:
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002) #this smoothing should get in right ballpark
        error = np.absolute(y_array - new_y_init)
        sm_error = df.gsmooth(x_array, error, var_y, .008)
        SNR = np.median(new_y_init / sm_error)

    vexp = 0.000205115 + 0.00699404*(SNR + 0.642253)**(-0.446162)

    if SNR > 150:
        vexp = .001

    return vexp, SNR

def deredden(a_v, r_v, wave, flux, var, model = 'f99'):
    wave = wave*u.Angstrom        # wavelengths
    flux = flux*u.Unit('W m-2 angstrom-1 sr-1')
    var = var*u.Unit('W^2 m-4 angstrom-2 sr-2')
    spec1d = Spectrum1D.from_array(wave, flux)
    spec1d_var = Spectrum1D.from_array(wave, var)
    red = ex.reddening(spec1d.wavelength, a_v = a_v, r_v = r_v, model=model)
    spec1d.flux *= red
    spec1d_var.flux *= (red**2.) #correct var too

    return spec1d.flux.value, spec1d_var.flux.value


def find_extrema(wavelength,sm_flux):
    #input a smoothed spectrum
    dw = 2.
    deriv_flux = np.diff(sm_flux)/dw
    abs_deriv = np.absolute(deriv_flux)
    extrema_locs = np.where(abs_deriv < .0001)
    maxima = []
    for loc in extrema_locs[0]:
        if loc != 0 and loc != len(deriv_flux)-1:
            if deriv_flux[loc-1] > 0 and deriv_flux[loc+1] < 0:
                maxima.append(loc)

    return maxima


def measure_si_ratio(wavelength,flux, varflux = None, vexp=.002, smooth=True, dm15 = None, plot=True):
    # if varflux == None:
    #   varflux = np.zeros(len(wavelength), float)
    if smooth:
        sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
    else:
        sm_flux = flux

    maxima = find_extrema(wavelength,sm_flux)
    max_waves = wavelength[maxima]
    max_fluxes = sm_flux[maxima]


    si_range_1 = np.where((wavelength > 5500.) & (wavelength < 5700.))
    si_range_2 = np.where((wavelength > 5700.) & (wavelength < 6000.))
    si_range_3 = np.where((wavelength > 6000.) & (wavelength < 6500.))

    # si_wave_1 = wavelength[si_range_1]
    # si_wave_2 = wavelength[si_range_2]
    # si_wave_3 = wavelength[si_range_3]

    # si_flux_1 = sm_flux[si_range_1]
    # si_flux_2 = sm_flux[si_range_2]
    # si_flux_3 = sm_flux[si_range_3]

    # si_max_1 = np.amax(si_flux_1)
    # si_max_2 = np.amax(si_flux_2)
    # si_max_3 = np.amax(si_flux_3)

    # si_max_index_1 = np.where(si_flux_1 == si_max_1)
    # si_max_index_2 = np.where(si_flux_2 == si_max_2)  
    # si_max_index_3 = np.where(si_flux_3 == si_max_3)

    # si_max_wave_1 = si_wave_1[si_max_index_1][0]
    # si_max_wave_2 = si_wave_2[si_max_index_2][0]
    # si_max_wave_3 = si_wave_3[si_max_index_3][0]

    m1s = []
    m2s = []
    m3s = []
    for m in maxima:
        if m in si_range_1[0]:
            m1s.append(m)
        if m in si_range_2[0]:
            m2s.append(m)
        if m in si_range_3[0]:
            m3s.append(m)

    if len(m1s) == 0 or len(m2s) == 0 or len(m3s) == 0:
        m1s = []
        m2s = []
        m3s = []
        sm_flux = df.gsmooth(wavelength, flux, varflux, vexp= vexp+.001)
        maxima = find_extrema(wavelength,sm_flux)
        for m in maxima:
            if m in si_range_1[0]:
                m1s.append(m)
            if m in si_range_2[0]:
                m2s.append(m)
            if m in si_range_3[0]:
                m3s.append(m)
    if len(m1s) == 0 or len(m2s) == 0 or len(m3s) == 0:
        print "Could not find maximum in a specified range!"
        return np.nan

    f1s = sm_flux[m1s]
    f2s = sm_flux[m2s]
    f3s = sm_flux[m3s]
    m1 = m1s[-1]
    m2 = m2s[-1]
    m3 = m3s[np.argmax(f3s)]

    # m1 = m1s[0]
    # m2 = m2s[-1]
    # m3 = m3s[-1]
    if dm15 != None and dm15 > 1.7:
        m1 = m1s[0]
        m2 = m2s[np.argmax(f2s)]

    si_max_wave_1 = wavelength[m1]
    si_max_wave_2 = wavelength[m2]
    si_max_wave_3 = wavelength[m3]

    si_max_1 = sm_flux[m1]
    si_max_2 = sm_flux[m2]
    si_max_3 = sm_flux[m3]

    weak_si_trough = np.where((wavelength >= si_max_wave_1) & (wavelength < si_max_wave_2))[0]
    strong_si_trough = np.where((wavelength >= si_max_wave_2) & (wavelength < si_max_wave_3))[0]

    interp_flux = copy.copy(flux)
    interp_flux[weak_si_trough] = np.interp(wavelength[weak_si_trough], [si_max_wave_1, si_max_wave_2], [si_max_1, si_max_2])
    interp_flux[strong_si_trough] = np.interp(wavelength[strong_si_trough], [si_max_wave_2, si_max_wave_3], [si_max_2, si_max_3])


    # v_strong, si_min_wave = measure_velocity(wavelength,flux, 5900., 6300.)
    # v_weak, si_weak_min_wave = measure_weak_si_velocity(wavelength,flux)

    # si_min_index = np.where(wavelength == si_min_wave)
    # si_weak_min_index = np.where(wavelength == si_weak_min_wave)

    #find max of diffs instead of finding minimum
    strong_line_diffs = interp_flux[strong_si_trough] - sm_flux[strong_si_trough]
    weak_line_diffs = interp_flux[weak_si_trough] - sm_flux[weak_si_trough]


    #Line ratio with fractional depths
    # strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])/interp_flux[si_min_index]
    # weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])/interp_flux[si_weak_min_index]

    # strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])
    # weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])

    strong_line = np.amax(strong_line_diffs)
    weak_line = np.amax(weak_line_diffs)
    strong_ind = np.where(strong_line_diffs == strong_line)
    weak_ind = np.where(weak_line_diffs == weak_line)

    ratio = weak_line/strong_line
    # ratio = ratio[0]


    # if ratio > .5:
    if plot:
        plt.plot(wavelength,flux)
        plt.plot(wavelength,sm_flux)
        plt.plot(wavelength, interp_flux)
        plt.plot(max_waves, max_fluxes, 'o', color='cyan')
        plt.plot(wavelength[strong_si_trough][strong_ind], interp_flux[strong_si_trough][strong_ind], 'o', color='orange')
        plt.plot(wavelength[strong_si_trough][strong_ind], sm_flux[strong_si_trough][strong_ind], 'o', color='orange')
        plt.plot(wavelength[weak_si_trough][weak_ind], interp_flux[weak_si_trough][weak_ind], 'o', color='orange')
        plt.plot(wavelength[weak_si_trough][weak_ind], sm_flux[weak_si_trough][weak_ind], 'o', color='orange')
        plt.xlim([5000.,7000.])
        plt.ylim([0.,.6])
        plt.show()

    # plt.plot(wavelength,flux)
    # plt.plot(wavelength,sm_flux)
    # plt.plot(wavelength, interp_flux)
    # plt.plot(wavelength[si_min_index], interp_flux[si_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_min_index], sm_flux[si_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_weak_min_index], interp_flux[si_weak_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_weak_min_index], sm_flux[si_weak_min_index], 'o', color='orange')
    # plt.xlim([5000.,7000.])
    # plt.show()

    return ratio

def measure_ca_ratio(wavelength,flux, varflux = None, wave1 = 3550., wave2=3680.):
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)

    # ca_range_1 = np.where((wavelength > 3550.) & (wavelength < 3680.))
    ca_range_1 = np.where((wavelength > wave1) & (wavelength < wave2))
    ca_range_2 = np.where((wavelength > 3890.) & (wavelength < 3970.))

    ca_wave_1 = wavelength[ca_range_1]
    ca_wave_2 = wavelength[ca_range_2]

    ca_flux_1 = sm_flux[ca_range_1]
    ca_flux_2 = sm_flux[ca_range_2]

    ca_max_1 = np.amax(ca_flux_1)
    ca_max_2 = np.amax(ca_flux_2)

    ca_max_index_1 = np.where(sm_flux == ca_max_1)
    ca_max_index_2 = np.where(sm_flux == ca_max_2)

    ca_max_wave_1 = wavelength[ca_max_index_1][0]
    ca_max_wave_2 = wavelength[ca_max_index_2][0]

    print ca_max_2/ca_max_1

    plt.plot(wavelength,flux)
    plt.plot(wavelength,sm_flux)
    plt.plot(wavelength[ca_max_index_1], sm_flux[ca_max_index_1], 'o', color='orange')
    plt.plot(wavelength[ca_max_index_2], sm_flux[ca_max_index_2], 'o', color='orange')
    plt.show()

    return ca_max_2/ca_max_1

def measure_verror(wavelength, flux, var_flux, wave1, wave2, n=100):
    vdist = []
    for i in range(n):
        sample_vexp = np.random.uniform(.001, .0045)
        sample_v, sample_si_min_wave, err = measure_velocity(wavelength, flux, wave1, wave2, clip=False, vexp=sample_vexp, plot=False, error=False)
        vdist.append(sample_v)

    sigma = np.std(vdist)

    return sigma

def calculate_velocity(wave, rest_wave):
    c = 299792.458
    v = -1.*c*((rest_wave/wave)**2. - 1)/(1+((rest_wave/wave)**2.))
    return v

def calculate_wave_from_velocity(velocity, rest_wave):
    c = 299792.458
    wave = rest_wave*np.sqrt((c+velocity)/(c-velocity))
    return wave

def measure_velocity(wavelength, flux, wave1, wave2, vexp=.001, clip=True, rest_wave=6355., varflux=None, plot=False, error=False, sn_name = None):

    sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)

    # sigclip_region = np.where((np.absolute(flux-sm_flux)> 3.*np.sqrt(varflux)))[0]
    old_wave = copy.deepcopy(wavelength)
    old_flux = copy.deepcopy(flux)
    # if clip:
    #     wavelength = np.delete(wavelength, sigclip_region)
    #     flux = np.delete(flux, sigclip_region)
    #     varflux = np.delete(varflux, sigclip_region)
    #     sm_flux = np.delete(sm_flux, sigclip_region)
    si_range = np.where((wavelength > wave1) & (wavelength < wave2))
    si_wave = wavelength[si_range]
    if len(si_wave) == 0:
        return np.nan, np.nan
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)

    # if len(si_min_index[0]) > 0. and (wavelength[-1] > wave2):
    if len(si_min_index[0]) > 0.:
        si_min_wave = si_wave[si_min_index][0]

        c = 299792.458 # km/s
        # rest_wave = 6355. #Angstroms

        v = c*((rest_wave/si_min_wave)**2. - 1)/(1+((rest_wave/si_min_wave)**2.))
        if error:
            sigma = measure_verror(wavelength, flux, varflux, wave1, wave2)

        if plot:
            f, ax1 = plt.subplots(1, 1, figsize=[7,7])
            zoom = (wavelength>wave1) & (wavelength < wave2)
            # zoom_old_sig = (old_wave[sigclip_region]>5500) & (old_wave[sigclip_region] < 6500)
            zoom_old = (old_wave>wave1) & (old_wave < wave2)
            norm = 1./np.amax(flux)
            plt.plot(old_wave[zoom_old], norm*old_flux[zoom_old], color='red')
            plt.plot(wavelength[zoom], norm*flux[zoom])
            plt.plot(wavelength[zoom], norm*sm_flux[zoom])
            if error:
                plt.plot(si_min_wave, norm*si_min, 'o', color='red', label = 'v = '+ str(np.round((-1.*v)/1000.,2)) + '+/- ' + str(np.round(sigma,2))+ ' 10^3 km/s')
            else:
                plt.plot(si_min_wave, norm*si_min, 'o', color='red', label = 'v = '+ str(np.round((-1.*v)/1000.,2)) + ' 10^3 km/s')
            plt.title(sn_name)
            plt.legend()
            # plt.xlim([5500.,6500.])
            # plt.ylim([np.median(flux[zoom])-.2,np.median(flux[zoom])+.2])
            # plt.savefig('../../Foundation/Spectra/new_new_new_spec/'+sn_name + '_' + str(np.round((-1.*v)/1000.,2))+'.png')
            plt.show()
    else:
        plt.plot(wavelength, flux, 'r')
        plt.show()
        v = np.nan
        si_min_wave = np.nan
        sigma = np.nan

    if error:
        return (-1.*v)/1000., si_min_wave, sigma
    else:
        return (-1.*v)/1000., si_min_wave, np.nan

def measure_vels(comps, sn_arrs, attr_name, boot_arrs = None, plot=False):
    avg_vs = []
    for comp in comps:
        v, si_min_wave = measure_velocity(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], 5800, 6400, plot=plot)
        avg_vs.append(v)
        print v

    if boot_arrs:
        errors = []
        boot_vel_arrs = []
        
        for boots in boot_arrs:
            b_vs = []
            for b in boots:
                v, si_min_wave = measure_velocity(b.wavelength[b.x1:b.x2], b.flux[b.x1:b.x2], 5800, 6400, plot=False)
                # print v
                if ~np.isnan(v) and v != None:
                    b_vs.append(v)
            boot_vel_arrs.append(b_vs)
            # errors.append(np.std(b_vs))
            plt.hist(b_vs)
            plt.show()
            
        low_errors = []
        up_errors = []
        for vlist in boot_vel_arrs:
            p = np.percentile(vlist, [18, 50, 82])
            vlower = p[1] - p[0]
            vupper = p[2] - p[1]
            low_errors.append(vupper)
            up_errors.append(vlower)
            print 'ERR: ', (vupper+vlower)/2

        errors = [low_errors, up_errors]
        
    avg_attrs = []
    for arr in sn_arrs:
        attr_list = []
        for i, spec in enumerate(arr):
            attr_list.append(spec.other_meta_data[attr_name])
        avg_attrs.append(np.average(attr_list))
    
    #Caveat: this includes combined spectra
    vels = []
    attrs = []
    attr_errs = []
    for arr in sn_arrs:
        for i, spec in enumerate(arr):
            v, si_min_wave = measure_velocity(spec.wavelength[spec.x1:spec.x2], spec.flux[spec.x1:spec.x2], 5800, 6300, plot=plot)
            attrs.append(spec.other_meta_data[attr_name])
            attr_errs.append(0)
            vels.append(v)
            print i,spec.name, v
    
    # if boot_arrs:
    #     plt.errorbar(avg_attrs, avg_vs, yerr=errors, fmt='o')
    # else:
    #     plt.scatter(avg_attrs, avg_vs)
    # plt.scatter(attrs, vels, color='orange')
    # plt.show()

    return boot_vel_arrs, avg_vs, avg_attrs, attrs, attr_errs, vels, errors

def plot_vels(vel_data, vwidth = .1, savename=None):
    boot_vel_arrs, avg_vs, avg_attrs, attrs, attr_errs, vels, errors = vel_data[0], vel_data[1], vel_data[2], vel_data[3], vel_data[4], vel_data[5], vel_data[6]
    plt.rc('font', family='serif')
    plt.figure(num = 1, dpi = 100, figsize = [7,7])
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    #     ax.get_yaxis().set_ticks([])
    #     plt.ylim([-311,12])
    plt.ylabel('Velocity ($10^3$ km s$^{-1}$)', fontsize = 20)
    # plt.xlabel('Hubble Residual (mag)', fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=20)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=10)

    plt.errorbar(attrs, vels, fmt='o', markersize=5, color = 'black', capsize=5, zorder=-10, label='Individual SNe')
    violin_parts = plt.violinplot(boot_vel_arrs, avg_attrs, points=40, widths=vwidth,
                      showmeans=False, showextrema=False, showmedians=False,
                      bw_method='silverman')
    
    print avg_vs
    diff = avg_vs[1] - avg_vs[0]
    # err = np.sqrt(errors[0][0]**2. + errors[1][1]**2.)
    # print 'Vdiff: ', diff, 'Err:', err, 'Sig: ', diff/err
    plt.errorbar(avg_attrs, avg_vs, fmt='*', markersize=20, color = 'orange', alpha=1., zorder=1, label= 'Composite \n Spectra')
#         vp = violin_parts['cmeans']
#         vp.set_edgecolor('darkblue')
#         vp.set_linewidth(4)
#         vp.set_alpha(1.)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('deepskyblue')
        pc.set_edgecolor('steelblue')
#             pc.set_color('red')
        pc.set_alpha(.5)
    # plt.xlim(-20,-7)
    # plt.axvline(x=0., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)
    # plt.text(-.45, -14.5, '+4 Days', fontsize=20)
    # plt.legend(bbox_to_anchor=(0.48, 0.45, 0.48, 0.5), fontsize=15)
    plt.gca().invert_yaxis()
    if savename is not None:
        plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()
#adapted from Rodrigo 
def max_wave(wavelength, flux, w1, w2, w3, vexp=.001, sys_error=False):
    
    sm_flux= df.gsmooth(wavelength, flux, None, vexp)

    wave_domain_1 = (wavelength > w1) & (wavelength < w2)
    elem_flux_1 = np.argmax(sm_flux[wave_domain_1]) #find minimum value within these flux vales to locate "dip
    max_wave_1 = wavelength[wave_domain_1][elem_flux_1] #find the corresponding wavelength
    
    wave_domain_2 = (wavelength > w2) & (wavelength < w3)
    elem_flux_2 = np.argmax(sm_flux[wave_domain_2]) #find minimum value within these flux vales to locate "dip
    max_wave_2 = wavelength[wave_domain_2][elem_flux_2] #find the corresponding wavelength

    if sys_error:
        #change continuum locations by a random even integer (at most 100 Angstroms)
        max_wave_1 = max_wave_1 + random.randint(-50, 50)*2.
        max_wave_2 = max_wave_2 + random.randint(-50, 50)*2.

    return max_wave_1, max_wave_2, wave_domain_1, wave_domain_2, elem_flux_1, elem_flux_2, sm_flux

def measure_EWs(sn_array, w1=7600., w2=8200., w3=9000., error=False):
    EWs = []
    phases = []
    errs = []
    morphs = []
    SNs = []
    for i, SN in enumerate(sn_array):
        if SN.wavelength[SN.x2] >= w3:
            vexp_auto, SNR = autosmooth(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], var_y=None)
            ew, domain = measure_EW(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], w1, w2, w3, vexp=vexp_auto, plot=True)
            if error:
                stat_err = ew_stat_error(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2], w1, w2, w3, domain, vexp=vexp_auto, num=50)
                sys_err = ew_sys_error(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2], w1, w2, w3, domain, vexp=vexp_auto, num=50)
            else:
                stat_err = np.nan
                sys_err = np.nan
            print i, ew, sys_err, stat_err, SN.phase
            # err_tot = np.sqrt(sys_err**2. + stat_err**2.)
            err_tot = np.nan # FIX THIS
            if not np.isnan(ew) and ew < 500. and stat_err < 100.:
                EWs.append(ew)
                phases.append(SN.phase)
                errs.append(err_tot)
                morphs.append(SN.ned_host)
                SN.nir_EW = ew
                SN.nir_EW_err = err_tot
                SNs.append(SN)
    mean_EW = np.nanmean(EWs)
    var_EW = np.nanstd(EWs)**2

    return mean_EW, var_EW, EWs, errs, phases, morphs, SNs

def measure_comp_diff_EW(boot_sn_arrays, w1=7600., w2=8200., w3=9000.):
    means = []
    varis = []
    EW_arrs = []
    for arr in boot_sn_arrays:
        mean_ew, var_EW, EWs = measure_EWs(arr, w1=w1, w2=w2, w3=w3)
        means.append(mean_ew)
        varis.append(var_EW)
        EW_arrs.append(EWs)

    diff = (means[0] - means[1])
    err = np.sqrt(varis[0]+varis[1])
    print 'Diff = ', diff, '+/-', err
    plt.hist(EW_arrs[0], color = 'blue')
    plt.hist(EW_arrs[1], color = 'red', alpha=.5)
    plt.show()

    return diff, err, means, varis, EWs

#adapted from Rodrigo 
def measure_EW(wavelength, flux, w1, w2, w3, vexp=.001, plot=False, sys_error=False):
    
    max_wave_1, max_wave_2, wave_domain_1, wave_domain_2, elem_flux_1, elem_flux_2, sm_flux = max_wave(wavelength, flux, w1, w2, w3, vexp=vexp, sys_error=sys_error)
    

    domain = (wavelength >= max_wave_1) & (wavelength <= max_wave_2)
    roi = (wavelength > w1) & (wavelength < w3)
    wave_range = wavelength[domain]
    flux_range = flux[domain]


    # line_elem = np.polyfit([max_wave_1, max_wave_2], [sm_flux[np.where(wavelength == max_wave_1)],
    #                                         sm_flux[np.where(wavelength == max_wave_2)]], 1)
    line_elem = np.polyfit([max_wave_1, max_wave_2], [sm_flux[wave_domain_1][elem_flux_1],
                                            sm_flux[wave_domain_2][elem_flux_2]], 1)

    line = line_elem[0] * wave_range + line_elem[1]
    norm = flux_range / line

    a_curve = np.trapz((norm), x = wave_range)
    a_line = max(wave_range) - min(wave_range)
    eq_width = a_line - a_curve
    
    if plot==True:
        plt.figure()
        plt.plot(wavelength[np.where((wavelength>w1) & (wavelength<w3))],
            flux[np.where((wavelength>w1) & (wavelength<w3))])
        plt.plot(wavelength[np.where((wavelength>w1) & (wavelength<w3))],
            sm_flux[np.where((wavelength>w1) & (wavelength<w3))], color='g')
#         plt.xlim([6200.,6400.])
        # plt.ylim([0.,3.])
        plt.axvline (x=max_wave_1, color = 'red')
        plt.axvline (x=max_wave_2, color = 'red')
        plt.plot(wavelength[domain],line, color='orange')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.show()
    
    return eq_width, roi

def ew_stat_error(wavelength, flux, ivar, w1, w2, w3, roi, vexp=.001, num=100):
    sm_flux = df.gsmooth(wavelength[roi], flux[roi], 1/ivar[roi], vexp)
    err = np.sqrt((1./ivar[roi]))
    # err = np.absolute(flux[roi] - sm_flux)
    sig = np.median(err)
    ews = []
    for i in range(0, num):
        new_flux = copy.deepcopy(sm_flux)
        delta_flux = np.random.normal(loc=0, scale=sig, size=len(new_flux))
        new_flux = new_flux + delta_flux
        ew, roi_ignore = measure_EW(wavelength[roi], new_flux, w1, w2, w3, vexp=vexp, plot=False)
        ews.append(ew)
    stat_err = np.std(ews)

    return stat_err

def ew_sys_error(wavelength, flux, ivar, w1, w2, w3, roi, vexp=.001, num=100):
    ews = []
    for i in range(0, num):
        vexp = random.uniform(0.001, 0.0045)
        ew, roi_ignore = measure_EW(wavelength[roi], flux[roi], w1, w2, w3, vexp=vexp, plot=False, sys_error=True)
        ews.append(ew)
    sys_err = np.std(ews)

    return sys_err


def measure_weak_si_velocity(wavelength, flux, varflux=None, plot=True):
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    varflux = np.zeros(len(wavelength), float)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    si_range = np.where((wavelength > 5560.) & (wavelength < 5820.))
    si_wave = wavelength[si_range]
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)
    si_min_wave = si_wave[si_min_index][0]

    c = 299792. # km/s
    si_rest_wave = 5979. #Angstroms

    v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
    # if plot:
    #   plt.plot(wavelength, flux)
    #   plt.plot(wavelength, sm_flux)
    #   plt.plot(si_min_wave, si_min, 'o')
    #   plt.xlim([5000.,7000.])
    #   plt.ylim([0.,.6])
    #   plt.show()
    return v, si_min_wave

    #carbon II is 6580A (2006bt good example)

def measure_si_velocity_from_raw(wavelength,flux,z, varflux=None):
    # varflux = 1./SN.ivar[SN.x1:SN.x2]
    # plt.plot(wavelength,flux)
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    wavelength = wavelength/(1.+z)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    si_range = np.where((wavelength > 5830.) & (wavelength < 6230.))
    si_wave = wavelength[si_range]
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)
    si_min_wave = si_wave[si_min_index][0]

    c = 299792.458 # km/s
    si_rest_wave = 6355. #Angstroms
    v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
    # v = c*(si_min_wave - si_rest_wave)/si_rest_wave
    # plt.plot(wavelength, flux)
    # plt.plot(wavelength, sm_flux)
    # plt.plot(si_min_wave, si_min, 'o', color='orange')
    # plt.show()
    return v, si_min_wave

def measure_C_velocity_from_raw(wavelength,flux,z, varflux=None):
    # varflux = 1./SN.ivar[SN.x1:SN.x2]
    # plt.plot(wavelength,flux)
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)
    wavelength = wavelength/(1.+z)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    C_range = np.where((wavelength > 6200.) & (wavelength < 6320.))
    C_wave = wavelength[C_range]
    C_flux = sm_flux[C_range]
    C_min = np.amin(C_flux)
    C_min_index = np.where(C_flux == C_min)
    C_min_wave = C_wave[C_min_index][0]

    c = 299792. # km/s
    C_rest_wave = 6580. #Angstroms
    v = c*((C_rest_wave/C_min_wave)**2. - 1)/(1+((C_rest_wave/C_min_wave)**2.))
    # v = c*(si_min_wave - si_rest_wave)/si_rest_wave
    # plt.plot(wavelength, flux)
    # plt.plot(wavelength, sm_flux)
    # plt.plot(C_min_wave, C_min, 'o', color='orange')
    # plt.show()
    return v, C_min_wave

def resample_spectrum(SN):
    wavelength = SN.wavelength[SN.x1:SN.x2]
    flux = SN.flux[SN.x1:SN.x2]
    varflux = 1./SN.ivar[SN.x1:SN.x2]


def measure_mags(wave, flux, filts = ['GROUND_JOHNSON_B','GROUND_JOHNSON_V']):
    lib = pyphot.get_library()
    mags = []
    for f in filts:
        func = lib[f]
        fflux = func.get_flux(wave, flux, axis = -1)
        mag = -2.5 * np.log10(fflux) - func.Vega_zero_mag
        mags.append(mag)
    return mags

def measure_comp_1m2(comps, filts = ['GROUND_JOHNSON_B','GROUND_JOHNSON_V'], boot_arrs = None, error=False):
    lib = pyphot.get_library()
    func1 = lib[filts[0]]
    func2 = lib[filts[1]]
    
    comp_1 = []
    comp_2 = []
    phases = []
    errors = []
    kpora.set_min_num_spec(comps, 5)
    for comp in comps:
        mag1, mag2 = measure_mags(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], filts = filts)
        comp_1.append(mag1)
        comp_2.append(mag2)
        phases.append(np.average(comp.phase_array[comp.x1:comp.x2]))
        
    if error:
        boot_1m2s = []
        errors = []
        for boots in boot_arrs:
            boot_1 = []
            boot_2 = []
            for b in boots:
                mag1, mag2 = measure_mags(b.wavelength[b.x1:b.x2], b.flux[b.x1:b.x2], filts=filts)
                if ~np.isnan(mag1) and ~np.isnan(mag2):
                    boot_1.append(mag1)
                    boot_2.append(mag2)
            diff = np.asarray(boot_1) - np.asarray(boot_2)
            boot_1m2s.append(diff)
        low_errors = []
        up_errors = []
        for clist in boot_1m2s:
            p = np.percentile(clist, [18, 50, 82])
            clower = p[1] - p[0]
            cupper = p[2] - p[1]
            low_errors.append(cupper)
            up_errors.append(clower)
        errors = [low_errors, up_errors]
        
    return phases, comp_1, comp_2, errors

def spectres(new_wavs, spec_wavs, spec_fluxes, spec_errs=None, fill=None,
             verbose=True):

    """
    Function for resampling spectra (and optionally associated
    uncertainties) onto a new wavelength basis.
    Parameters
    ----------
    new_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired for the
        spectrum or spectra.
    spec_wavs : numpy.ndarray
        1D array containing the current wavelength sampling of the
        spectrum or spectra.
    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths specified in
        spec_wavs, last dimension must correspond to the shape of
        spec_wavs. Extra dimensions before this may be used to include
        multiple spectra.
    spec_errs : numpy.ndarray (optional)
        Array of the same shape as spec_fluxes containing uncertainties
        associated with each spectral flux value.
    fill : float (optional)
        Where new_wavs extends outside the wavelength range in spec_wavs
        this value will be used as a filler in new_fluxes and new_errs.
    verbose : bool (optional)
        Setting verbose to False will suppress the default warning about
        new_wavs extending outside spec_wavs and "fill" being used.
    Returns
    -------
    new_fluxes : numpy.ndarray
        Array of resampled flux values, first dimension is the same
        length as new_wavs, other dimensions are the same as
        spec_fluxes.
    new_errs : numpy.ndarray
        Array of uncertainties associated with fluxes in new_fluxes.
        Only returned if spec_errs was specified.
    """

    # Rename the input variables for clarity within the function.
    old_wavs = spec_wavs
    old_fluxes = spec_fluxes
    old_errs = spec_errs

    # Arrays of left hand sides and widths for the old and new bins
    old_lhs = np.zeros(old_wavs.shape[0])
    old_widths = np.zeros(old_wavs.shape[0])
    old_lhs = np.zeros(old_wavs.shape[0])
    old_lhs[0] = old_wavs[0]
    old_lhs[0] -= (old_wavs[1] - old_wavs[0])/2
    old_widths[-1] = (old_wavs[-1] - old_wavs[-2])
    old_lhs[1:] = (old_wavs[1:] + old_wavs[:-1])/2
    old_widths[:-1] = old_lhs[1:] - old_lhs[:-1]
    old_max_wav = old_lhs[-1] + old_widths[-1]

    new_lhs = np.zeros(new_wavs.shape[0]+1)
    new_widths = np.zeros(new_wavs.shape[0])
    new_lhs[0] = new_wavs[0]
    new_lhs[0] -= (new_wavs[1] - new_wavs[0])/2
    new_widths[-1] = (new_wavs[-1] - new_wavs[-2])
    new_lhs[-1] = new_wavs[-1]
    new_lhs[-1] += (new_wavs[-1] - new_wavs[-2])/2
    new_lhs[1:-1] = (new_wavs[1:] + new_wavs[:-1])/2
    new_widths[:-1] = new_lhs[1:-1] - new_lhs[:-2]

    # Generate output arrays to be populated
    new_fluxes = np.zeros(old_fluxes[..., 0].shape + new_wavs.shape)

    if old_errs is not None:
        if old_errs.shape != old_fluxes.shape:
            raise ValueError("If specified, spec_errs must be the same shape "
                             "as spec_fluxes.")
        else:
            new_errs = np.copy(new_fluxes)

    start = 0
    stop = 0

    # Calculate new flux and uncertainty values, looping over new bins
    for j in range(new_wavs.shape[0]):

        # Add filler values if new_wavs extends outside of spec_wavs
        if (new_lhs[j] < old_lhs[0]) or (new_lhs[j+1] > old_max_wav):
            new_fluxes[..., j] = fill

            if spec_errs is not None:
                new_errs[..., j] = fill

            if (j == 0) and verbose:
                print("\nSpectres: new_wavs contains values outside the range "
                      "in spec_wavs. New_fluxes and new_errs will be filled "
                      "with the value set in the 'fill' keyword argument (nan "
                      "by default).\n")
            continue

        # Find first old bin which is partially covered by the new bin
        while old_lhs[start+1] <= new_lhs[j]:
            start += 1

        # Find last old bin which is partially covered by the new bin
        while old_lhs[stop+1] < new_lhs[j+1]:
            stop += 1

        # If new bin is fully inside an old bin start and stop are equal
        if stop == start:
            new_fluxes[..., j] = old_fluxes[..., start]
            if old_errs is not None:
                new_errs[..., j] = old_errs[..., start]

        # Otherwise multiply the first and last old bin widths by P_ij
        else:
            start_factor = ((old_lhs[start+1] - new_lhs[j])
                            / (old_lhs[start+1] - old_lhs[start]))

            end_factor = ((new_lhs[j+1] - old_lhs[stop])
                          / (old_lhs[stop+1] - old_lhs[stop]))

            old_widths[start] *= start_factor
            old_widths[stop] *= end_factor

            # Populate new_fluxes spectrum and uncertainty arrays
            f_widths = old_widths[start:stop+1]*old_fluxes[..., start:stop+1]
            new_fluxes[..., j] = np.sum(f_widths, axis=-1)
            new_fluxes[..., j] /= np.sum(old_widths[start:stop+1])

            if old_errs is not None:
                e_wid = old_widths[start:stop+1]*old_errs[..., start:stop+1]

                new_errs[..., j] = np.sqrt(np.sum(e_wid**2, axis=-1))
                new_errs[..., j] /= np.sum(old_widths[start:stop+1])

            # Put back the old bin widths to their initial values
            old_widths[start] /= start_factor
            old_widths[stop] /= end_factor

    # If errors were supplied return both new_fluxes and new_errs.
    if old_errs is not None:
        return np.array([new_wavs, new_fluxes, new_errs])
        # return new_fluxes, new_errs

    # Otherwise just return the new_fluxes spectrum array
    else:
        # return new_fluxes
        return np.array([new_wavs, new_fluxes])

def process_keck_file(spec):
    meta_dict = {}
    spec_fits=fits.open(spec)
    
    filename = spec.split('/')[-1]
    head=spec_fits[0].header
    meta_dict['FILENAME']  = filename
    meta_dict['OBJECT']  = head['OBJECT']
    
    flux, err = np.transpose(spec_fits[0].data)

    snr = np.round(np.median(flux/err),3)
    meta_dict['SNR'] = snr
    
    wavezero=float(head['CRVAL1'])
    wavedelt=float(head['CDELT1'])
    wave=np.arange(len(flux))*wavedelt+wavezero
    airmass=float(head['AIRMASS'])
    
    meta_dict['EXPTIME']  = float(head['EXPTIME'])
    meta_dict['AIRMASS'] = np.round(float(head['AIRMASS']),3)
    meta_dict['MINWAVE'] = np.round(float(head['W_RANGE'].split()[0]),4)
    meta_dict['MAXWAVE'] = np.round(float(head['W_RANGE'].split()[1]),4)
    meta_dict['WAVEDELT'] = np.round(wavedelt,4)
    meta_dict['CRVAL'] = np.round(wavezero,4)
    meta_dict['MJD'] = float(head['MJD-OBS'])
    meta_dict['UTC']  = head['UTC']
    meta_dict['POS_ANG'] = np.round(float(head['ROTPOSN'])+90, 2)
    meta_dict['FLUX_OBJ'] = head['FLUX_OBJ']
    
    radec = SkyCoord(head['RA'], head['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
    meta_dict['RA_OBS'] = np.round(radec.ra.deg, 4)
    meta_dict['DEC_OBS'] = np.round(radec.dec.deg, 4)
#     apo = Observer.at_site("Keck")
#     date = head['DATE_BEG']
#     par_ang = apo.parallactic_angle(date.split('T')[0] + ' ' + date.split('T')[1], radec)
#     print (par_ang)
    
    meta_dict['OBSERVER'] = head['OBSERVER'].strip()
    meta_dict['REDUCER'] = head['REDUCER']
    
    meta_dict['INSTRUMENT'] = head['INSTRUME'].strip()
    meta_dict['GRATING'] = head['GRANAME'].strip()
    meta_dict['GRISM'] = head['GRISNAME'].strip()
    meta_dict['DICHROIC'] = head['DICHNAME'].strip()
    meta_dict['SLIT'] = head['SLITNAME'].strip()
    meta_dict['TELESCOPE'] = head['TELESCOP'].strip()
    
    if 'combined' in filename:
        meta_dict['COMBINED'] = 1
    else:
        meta_dict['COMBINED'] = 0
        
    if 'BAD' in filename:
        meta_dict['AP_GREATER_SEEING'] = 0
    else:
        meta_dict['AP_GREATER_SEEING'] = 1
        
    if 'arcsec' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('arcsec')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        if 'SN' in filename.split('arcsec')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'kpc' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('kpc')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'kpc'
        if 'SN' in filename.split('kpc')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'rkron' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('rkron')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 1
    else:
        meta_dict['AP_SIZE'] = None
        meta_dict['AP_UNIT'] = None
        meta_dict['AP_LOC'] = None
        meta_dict['IS_KRON_RAD'] = 0
        
    return wave, flux, err, meta_dict
    
def process_lick_file(spec):
    meta_dict = {}
    spec_fits=fits.open(spec)
    
    filename = spec.split('/')[-1]
    head=spec_fits[0].header
    meta_dict['FILENAME']  = filename
    meta_dict['OBJECT']  = head['OBJECT']
    
    wavezero=float(head['CRVAL1'])
    wavedelt=float(head['CD1_1'])
    wave=np.arange(len(spec_fits[0].data[0][0]))*wavedelt+wavezero

    norm = matplotlib.colors.Normalize(vmin=0,vmax=len(spec_fits[0].data[0]))
    c_m = matplotlib.cm.plasma
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    labels = ['AP: -5.5:5.5',
              'AP: -8:-6',
              'AP: -6:-4',
              'AP: -4:-2',
              'AP: -2:0',
              'AP: 0:2',
              'AP: 2:4',
              'AP: 4:6',
              'AP: 6:8',
              ]
    labels = ['none',
              'none',
              '-6:-4/-8:-6',
              '-4:-2/-6:-4',
              '-2:0/-4:-2',
              '0:2/-2:0',
              '2:4/0:2',
              '4:6/2:4',
              '6:8/4:6'
              ]
    # for i, flux in enumerate(spec_fits[0].data[0]):
    #     roi = (wave > 6000) & (wave < 7000)
    #     sca = np.median(spec_fits[0].data[0][0][roi])/np.median(flux[roi])
    #     plt.plot(wave,sca*flux, color=s_m.to_rgba(i), label = labels[i])
    # plt.legend()
    # plt.show()

    # for i, flux in enumerate(spec_fits[0].data[0]):
    #     roi = (wave > 6000) & (wave < 7000)
    #     sca = np.median(spec_fits[0].data[0][0][roi])/np.median(flux[roi])
    #     plt.plot(wave,sca*flux/spec_fits[0].data[0][0], color=s_m.to_rgba(i), label = labels[i])
    # plt.legend()
    # plt.show()

    for i, flux in enumerate(spec_fits[0].data[0]):
        if i > 1:
            roi = (wave > 5800) & (wave < 6000)
            sca = np.median(spec_fits[0].data[0][i-1][roi])/np.median(flux[roi])
            plt.plot(wave,sca*flux/spec_fits[0].data[0][i-1], color=s_m.to_rgba(i), label = labels[i])
    plt.legend()
    plt.show()
    raise TypeError

    flux, err = np.transpose(spec_fits[0].data)
    snr = np.round(np.median(flux/err),3)
    meta_dict['SNR'] = snr
    
    wavezero=float(head['CRVAL1'])
    wavedelt=float(head['CDELT1'])
    wave=np.arange(len(flux))*wavedelt+wavezero
    airmass=float(head['AIRMASS'])
    
    meta_dict['EXPTIME']  = float(head['EXPTIME'])
    meta_dict['AIRMASS'] = np.round(float(head['AIRMASS']),3)
    meta_dict['MINWAVE'] = np.round(float(head['W_RANGE'].split()[0]),4)
    meta_dict['MAXWAVE'] = np.round(float(head['W_RANGE'].split()[1]),4)
    meta_dict['WAVEDELT'] = np.round(wavedelt,4)
    meta_dict['CRVAL'] = np.round(wavezero,4)
    ut_date = head['DATE-OBS'].strip()
    t = Time(ut_date, format='isot')
    meta_dict['MJD'] = float(t.mjd)
    meta_dict['UTC']  = ut_date.split('T')[1]
    meta_dict['POS_ANG'] = np.round(float(head['TUB']), 2) #TODO:check that tub is correct
    meta_dict['FLUX_OBJ'] = head['FLUX_OBJ']
    
    radec = SkyCoord(head['RA'], head['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
    meta_dict['RA'] = np.round(radec.ra.deg, 4)
    meta_dict['DEC'] = np.round(radec.dec.deg, 4)
#     apo = Observer.at_site("Keck")
#     date = head['DATE_BEG']
#     par_ang = apo.parallactic_angle(date.split('T')[0] + ' ' + date.split('T')[1], radec)
#     print (par_ang)
    
    meta_dict['OBSERVER'] = head['OBSERVER'].strip()
    meta_dict['REDUCER'] = head['REDUCER']
    
    meta_dict['INSTRUMENT'] = head['VERSION'].strip()
    meta_dict['GRATING'] = head['GRATNG_N'].strip()
    meta_dict['GRISM'] = head['GRISM_N'].strip()
    meta_dict['DICHROIC'] = head['BSPLIT_N'].strip()
    meta_dict['SLIT'] = head['SLIT_N'].strip()
    meta_dict['TELESCOPE'] = head['OBSERVAT'].strip()
    
    if 'combined' in filename:
        meta_dict['COMBINED'] = 1
    else:
        meta_dict['COMBINED'] = 0
        
    if 'BAD' in filename:
        meta_dict['AP_GREATER_SEEING'] = 0
    else:
        meta_dict['AP_GREATER_SEEING'] = 1
        
    if 'arcsec' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('arcsec')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        if 'SN' in filename.split('arcsec')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'kpc' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('kpc')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'kpc'
        if 'SN' in filename.split('kpc')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'rkron' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('rkron')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 1
    else:
        meta_dict['AP_SIZE'] = np.nan
        meta_dict['AP_UNIT'] = 'None'
        meta_dict['AP_LOC'] = 'None'
        meta_dict['IS_KRON_RAD'] = 0
        
    return wave, flux, err, meta_dict

def process_soar_file(head):
    return

line_dict = {'H':       ([6562.79, 4861.35, 4340.472, 4101.734], 'mediumblue'),
             'He':      ([5876.], 'black'),
             '[O III]': ([4958.911, 5006.843, 4363.210], 'magenta'),
             '[O II]':  ([3726.032, 3728.815], 'magenta'),
             '[N II]':  ([6548.050, 6583.460], 'darkorange'),
             '[S II]':  ([6716.440, 6730.810], 'darkgreen'),
             'Ca H':    ([3968.5], 'red'),
             'Ca K':    ([3933.7], 'red'),
             'G-band':  ([4304.4], 'gold'),
             'Mg':      ([5175.3], 'purple'),
             'Na ID':   ([5894.0], 'lime'),
             '[Fe II]': ([7155.1742], 'slategray'),
             '[Ca II]': ([7291.47, 7323.89], 'indigo'),
             '[Ni II]': ([7377.83], 'coral')
             }

if __name__ == "__main__":

    description = "For investigation of individual spectra"
    usage = "%prog    \t [option] \n Recommended syntax: %prog"
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option("-v", "--velocity", dest="vel", action="store_true",
                      help='Measure a line velocity from an absorption minimum')
    parser.add_option("-n", "--plot-neb-lines", dest="neblines", action="store_true",
                      help='Plot nebular emission lines')
    parser.add_option("-l", "--plot-all-lines", dest="alllines", action="store_true",
                      help='Plot all emission lines')
    parser.add_option("-c", "--csv", dest="csv", action="store_true",
                      help='Spectrum is a csv file')
    parser.add_option("-s", "--scale", dest="scale", action="store_true",
                      help='Scale to peak in spectrum')
    parser.add_option("-i", "--interp", dest="interp", action="store_true",
                      help='Interpolate the spectrum')

    option, args = parser.parse_args()
    _vel= option.vel
    _neblines= option.neblines
    _lines= option.alllines
    _csv= option.csv
    _sca= option.scale
    _interp= option.interp

    c = 299792.458
    # v = c*((6355./6012.)**2. - 1)/(1+((6355./6012.)**2.))
    # v = c*((6562.79/6453.6)**2. - 1)/(1+((6562.79/6453.6)**2.))
    # v = c*((5876/5857.85)**2. - 1)/(1+((5876/5857.85)**2.))
    # print v

    spec_file_names = raw_input("Choose an fits/ascii/csv spectrum file: ")
    spec_files = []
    for spec_file in spec_file_names.split(' '):
        spec_files.append(spec_file)

    # spec_files = ['sn2021fxy-20210405_UCB.flm',
    #               '2021fxy-red-20210406_ap1_center.flm',
    #               '2021fxy-red-20210406_ap1_right.flm',
    #               '2021fxy-red-20210406_ap1_left.flm'
    #               ]

    # spec_files = ['d2021fxy_kast_red_1_ex.fits']
    # spec_files = ['dBD262606_kast_red_1_ex.fits']

    # binning comparison
    # spec_files = ['sn2006jc-20061123.649-br.flm',
    #               'sn2006jc-blue-20061123_ap1_newarc_spectres.flm',

    #             'sn2006jc-combined-20061123_ap1_newarc_spectres.flm'
    #             ]


    # spec_files = ['sn2006lv-20061123.657-br.flm',
    #               'sn2006lv-combined-20061123_ap1.flm']
    # spec_files = ['sn2006my-20061123.644-br.flm',
    #               'sn2006my-combined-20061123_ap1.flm',
    #               'sn2006my-combined-20061123_ap1_noflat.flm']

    # spec_files = ['sn2006my-20061123.644-br.flm',
    #               'sn2006my-blue-20061123_ap1.flm',
    #               'sn2006my-blue-20061123_ap1_diff_std.flm'
    #               # 'sn2006my-blue-20061123_ap1_newshift.flm'
    #               ]

    # best arc solution
    # spec_files = ['sn2006my-20061123.644-br.flm',
    #               'sn2006my-combined-20061123_ap1_spectres.flm',
    #               'sn2006my-combined-20061123_ap1_newarc.flm',
    #               'sn2006my-combined-20061123_ap1_new2.flm'
    #               ]

    # full wave range
    # spec_files = ['sn2006my-20061123.644-br.flm',
    #           'sn2006my-combined-20061123_ap1_ashrebin.flm',
    #           'sn2006my-combined-20061123_ap1_new_arc.flm'
    #           # 'sn2006my-blue-20061123_ap1_newshift.flm'
    #           ]
    # spec_files = ['sn2005hk-20061123.234-br.flm',
    #               'sn2005hk-combined-20061123_ap1.flm']
    # spec_files = ['sn2006or-20061123.654-br.flm',
    #               'sn2006or-combined-20061123_ap1.flm']
    # spec_files = ['sn2006x-20061123.638-br.flm',
    #               'sn2006XoldSNIa-combined-20061123_ap1.flm']


    # spec_files = ['r185-20061123.291-br.flm',
    #               'r185-combined-20061123_ap1.flm']


    # spec_files = ['r206-20061123.260-br.flm',
    #               'r206-combined-20061123_ap1.flm']
    # spec_files = ['r193-20061123.422-br.flm',
    #               'r193-combined-20061123_ap1_newarc.flm']
    # # spec_files = ['sn2006ce-20061123.357-br.flm',
    # #               'sn2006ce-combined-20061123_ap1.flm',
    # #               'sn2006ce-combined-20061123_ap1_w_hz44.flm']
    # spec_files = ['sn2006ce-20061123.357-br.flm',
    #               'sn2006ce-combined-20061123_ap1.flm'
    #               # 'sn2006ce-combined-20061123_ap1_noflat.flm'
    #               ]

    # spec_files = [
    #               'sn2006ce-20061123.357-br.flm',
    #               'sn2006ce-blue-20061123_ap1.flm',
    #               'sn2006ce-blue-20061123_ap1_noflat.flm']

    # spec_files = ['2017hmf-combined-20180611_ap6_3.0_kpc.fits',
    #               '2017hmf-combined-20180711_ap1.fits']
    # spec_files = ['SN2018bdm-combined-20180611_ap1_georgios.fits',
    #               '2018bdm-combined-20180711_ap1.fits']
    # msfiles = ['dBD262606_1_lris_red_1_ex0610.fits', 'dBD262606_lris_red_1_ex.fits']
    # msfiles = ['dMIRAPMFM-350_lris_red_1_ex.fits', 'dBD262606_07_lris_red_1_ex.fits']
    
    # waves = []
    # fluxes = []
    # for i, msfile in enumerate(msfiles):
    #     multifits=fits.open(msfile)
    #     multispec=multifits[0].data
    #     mshead=multifits[0].header

    #     crval=float(mshead['CRVAL1'])
    #     cdelt=float(mshead['CD1_1'])
    #     npix=float(mshead['NAXIS1'])
    #     wave=np.arange(npix)*cdelt + crval
    #     flux = multispec[1,0,:]
    #     sca = 1.
    #     if _sca:
    #         roi = (wave > 5800) & (wave < 6200)
    #         sca = 1./np.median(flux[roi])
    #     plt.plot(wave, sca*flux)
    #     waves.append(wave)
    #     fluxes.append(sca*flux)
    # plt.show()


    # msfiles = ['dBD174708_lris_red_1_ex.fits', 'dBD274708_1_lris_red_1_ex0610.fits']
    # waves_17 = []
    # fluxes_17 = []
    # for i, msfile in enumerate(msfiles):
    #     multifits=fits.open(msfile)
    #     multispec=multifits[0].data
    #     mshead=multifits[0].header

    #     crval=float(mshead['CRVAL1'])
    #     cdelt=float(mshead['CD1_1'])
    #     npix=float(mshead['NAXIS1'])
    #     wave=np.arange(npix)*cdelt + crval
    #     flux = multispec[1,0,:]
    #     sca = 1.
    #     if _sca:
    #         roi = (wave > 5800) & (wave < 6200)
    #         sca = 1./np.median(flux[roi])
    #     plt.plot(wave, sca*flux)
    #     waves_17.append(wave)
    #     fluxes_17.append(sca*flux)
    # plt.show()



    # flux1 = np.interp(waves[0], waves[1], fluxes[1])
    # plt.plot(waves[0], np.asarray(fluxes[0])/flux1, label='BD284211')

    # # flux1_17 = np.interp(waves_17[0], waves_17[1], fluxes_17[1])
    # # plt.plot(waves_17[0], np.asarray(fluxes_17[0])/flux1_17,label='BD174708')

    # plt.legend(loc=2)
    # plt.show()
    # raise TypeError

    # spec_files = ['2021bls-combined-20210211_ap1.flm',
    #               'SN2016hnk_osc.json']
    # spec_files = ['2021bls-combined-20210211_ap1.flm',
    #               '2020wnt-combined-20210211_ap1.flm',]
    # spec_labels = ['2021bls', '2016hnk +4 days']
    # spec_labels = ['2021bls', '2020wnt']
    waves = []
    fluxes = []
    for j, spec_file in enumerate(spec_files):
        if spec_file.endswith('.flm'):
            data = np.genfromtxt(spec_file, unpack=True)
        elif _csv or spec_file.endswith('.csv'):
            data = np.genfromtxt(spec_file, unpack=True, delimiter=',')
        # elif spec_file.endswith('.json'):
        #     with open(spec_file) as f:
        #         data_json = json.load(f)
        #         wave = np.transpose(data_json['SN2016hnk']['spectra'][0]['data'])[0].astype(np.float)
        #         flux = np.transpose(data_json['SN2016hnk']['spectra'][0]['data'])[1].astype(np.float)
        #         err = None
        #         data = [wave, flux]
        else:
            tscope = raw_input("Which telescope? [keck]: ") or 'keck'
            if tscope == 'keck':
                wave, flux, err, meta_dict = process_keck_file(spec_file)
            elif tscope == 'lick':
                wave, flux, err, meta_dict = process_lick_file(spec_file)
            data = [wave,flux]

        redshift = raw_input("Redshift [0]: ") or 0.
        redshift=float(redshift)
        wavelength = data[0]
        flux = data[1]

        if _interp:
            dw = raw_input("New A/pix? [2]: ") or 2.
            dw = float(dw)
            # interp_wave = np.arange(math.ceil(wavelength[0])+1.*dw, math.floor(wavelength[-1])-1.*dw, dtype=float, step=dw)
            interp_wave = np.arange(3200, 9100, dtype=float, step=dw)
            # interp_wave = np.arange(3200, 5600, dtype=float, step=dw)
            binned_data = spectres(interp_wave, wavelength, flux, spec_errs=None, fill=None, verbose=True)
            wavelength = binned_data[0]
            flux = binned_data[1]
            # np.savetxt('20esm_full_first_pass_binned.flm', np.transpose([wavelength,flux]))

        wavelength = wavelength/(1.+redshift)

        # data_12dn = np.genfromtxt('2012dn_m14.flm', delimiter=',', unpack=True)
        # # print data_12dn
        # redshift_12dn = 0.010187
        # wavelength_12dn = data_12dn[0]
        # flux_12dn = data_12dn[1]
        # wavelength_12dn = wavelength_12dn/(1.+redshift_12dn)

        # plt.plot(wavelength_12dn,flux_12dn*(10**15))
        sca = 1.
        if _sca:
            roi = (wavelength > 4000) & (wavelength < 5300)

            # roi = (wavelength > 6000) & (wavelength < 7000)
            if len(flux[roi]) == 0:
                roi = (wavelength > 6000) & (wavelength < 7000)
                roi_0 = (waves[0] > 6000) & (waves[0] < 7000)
                # print fluxes[0]
                sca = np.median(fluxes[0][roi_0])/np.median(flux[roi])
            else:   
                sca = 1./np.median(flux[roi])
            
            # sca = 1./np.median(flux[0:10])
        # plt.plot(wavelength,sca*flux, linewidth=1, drawstyle='steps-mid', label=spec_labels[j])
        plt.plot(wavelength,sca*flux, linewidth=1, drawstyle='steps-mid', label=spec_file)
        waves.append(wavelength)
        fluxes.append(sca*flux)


        if _vel:
            r_wave = raw_input("Rest Wavelength [6355]: ") or 6355.
            r_wave = float(r_wave)
            wave_range = raw_input("Wavelength Range [5800 6400]: ") or '5800 6400'
            wave1 = float(wave_range.split()[0])
            wave2 = float(wave_range.split()[1])
            vexp, SNR = find_vexp(wavelength, flux)
            vel_data = measure_velocity(wavelength, flux, wave1, wave2, vexp=vexp, clip=True, rest_wave=r_wave, varflux=None, plot=True, error=False)
            print vel_data

    # if _neblines:
    #     color = ['mediumblue', 'magenta', 'magenta', 'darkorange']
    #     labels = ['Fe II', 'Ca II', '', 'Ni II']
    #     line_waves = [7155.1742, 7291.47, 7323.89, 7377.83]

    #     for i, line in enumerate(line_waves):
    #         if i != 2:
    #             plt.axvline(line, color=color[i], linestyle = ':', linewidth=3, label=labels[i])
    #         else:
    #             plt.axvline(line, color=color[i], linestyle = ':', linewidth=3)
    # plt.legend()
    # plt.show()

    if _lines:

        for i, line_name in enumerate(line_dict):
            for line in line_dict[line_name][0]:
                plt.axvline(line, color=line_dict[line_name][1], linestyle = ':', linewidth=3)
                plt.text(line+5, np.amax(sca*fluxes[0]) - .05*np.amax(sca*fluxes[0]), line_name, fontsize = 15, horizontalalignment='center', rotation=90)
    plt.legend()
    plt.show()

    #for comparing spectra
    # for i,flux in enumerate(fluxes):
    #     plt.plot(wavelength,fluxes[i]/fluxes[0], linewidth=1, drawstyle='steps-mid')
    #     plt.axhline(1, color='k', linewidth=2)
    #     plt.axhline(1.05, color='r', linewidth=2)
    #     plt.axhline(.95, color='r', linewidth=2)
    #     plt.axhline(1.01, color='g', linewidth=2)
    #     plt.axhline(.99, color='g', linewidth=2)
    # plt.show()
    # for i,flux in enumerate(fluxes):
    #     plt.plot(wavelength,fluxes[i]-fluxes[0], linewidth=1, drawstyle='steps-mid')
    #     plt.axhline(0, color='k', linewidth=3)
    # plt.show()




