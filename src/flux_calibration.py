import matplotlib
import matplotlib.pyplot as plt
import kaepora as kpora
import kaepora_plot as kplot
import flux_calibration as fc
import pyphot
import numpy as np
from scipy.interpolate import splrep, splev
import scipy.optimize as opt
from specutils import extinction as ex
from specutils import Spectrum1D
from astropy import units as u
import sys,os
import copy

# color_dict = {'U': 'magenta',
#               'B': 'blue',
#               'V': 'green',
#               'R': 'red',
#               'I': 'brown',

#               'u': 'purple',
#               'g': 'seagreen',
#               'r': 'crimson',
#               'i': 'orangered',
#               'z': 'darkred',

#               "u'": 'purple',
#               "g'": 'seagreen',
#               "r'": 'crimson',
#               "i'": 'orangered',
#               "z'": 'darkred',

#               'H': 'yellow',
#               'K': 'teal',
#               'J': 'crimson',
#               'V': 'orange',
#               'Y': 'pink',

#               'W1': 'gray',
#               'W2': 'black',
#               'M2': 'darkgray',

#               'C': 'limegreen',
#               'Js': 'gold',
#               'Ks': 'lightgray'}

color_dict = {'U': 'magenta',
          'B': 'blue',
          'V': 'green',
          'R': 'red',
          'I': 'brown',

          'u': 'purple',
          'g': 'seagreen',
          'r': 'crimson',
          'i': 'orangered',
          'z': 'darkred',

          "u'": 'purple',
          "g'": 'seagreen',
          "r'": 'crimson',
          "i'": 'orangered',
          "z'": 'darkred',

          'H': 'yellow',
          'K': 'teal',
          'J': 'crimson',
          'V': 'orange',
          'V0': 'orange',
          'Y': 'pink',

          'W1': 'gray',
          'W2': 'black',
          'M2': 'darkgray',

          'C': 'limegreen',
          'Js': 'gold',
          'Ks': 'lightgray',
          'Jrc2': 'lavender',
          'Hdw': 'darkblue',
          'Ydw': 'darkgreen',
          'Jdw': 'darkgray',
             }

def get_band_data(phot, band):
    mjds = []
    mags = []
    errs = []

    telescope_dict = {"None": 0}
    system_dict = {"None": 0}
    observatory_dict = {"None": 0}
    instrument_dict = {"None": 0}

    for mag in phot[band][1]:

        if mag[1].get('telescope', None) is not None:
            tscope = mag[1].get('telescope', None)
            if tscope not in telescope_dict:
                telescope_dict[tscope] = 1
            else:
                telescope_dict[tscope] = telescope_dict[tscope] + 1
        else:
            telescope_dict['None'] = telescope_dict['None'] + 1

        if mag[1].get('system', None) is not None:
            sys = mag[1].get('system', None)
            if sys not in system_dict:
                system_dict[sys] = 1
            else:
                system_dict[sys] = system_dict[sys] + 1
        else:
            system_dict['None'] = system_dict['None'] + 1

        if mag[1].get('observatory', None) is not None:
            obs = mag[1].get('observatory', None)
            if obs not in observatory_dict:
                observatory_dict[obs] = 1
            else:
                observatory_dict[obs] = observatory_dict[obs] + 1
        else:
            observatory_dict['None'] = observatory_dict['None'] + 1

        if mag[1].get('instrument', None) is not None:
            inst = mag[1].get('instrument', None)
            if inst not in instrument_dict:
                instrument_dict[inst] = 1
            else:
                instrument_dict[inst] = instrument_dict[inst] + 1
        else:
            instrument_dict['None'] = instrument_dict['None'] + 1
    print band
    print telescope_dict
    print system_dict
    print observatory_dict
    print instrument_dict

    max_sys = max(system_dict, key=lambda k: system_dict[k])
    print max_sys

    for i, mag in enumerate(phot[band][1]):
        sys = mag[1].get('system', "None")
        if  sys == max_sys:
            mjds.append(float(phot[band][0][i]))
            mags.append(float(mag[0]))
            if 'e_magnitude' in mag[1]:
                errs.append(float(mag[1]['e_magnitude']))
            elif 'upperlimit' in mag[1] and mag[1]['upperlimit']:
                errs.append(np.nan)
            else:
                errs.append(0.5)
    print len(mjds)

    mjds = np.asarray(mjds)
    mags = np.asarray(mags)
    errs = np.asarray(errs)

    detection_inds = np.where(~np.isnan(errs))[0]
    mjds = mjds[detection_inds]
    mags = mags[detection_inds]
    errs = errs[detection_inds]
    print
    return mjds, mags, errs


def interp_LC(mjds, mags, errs, s = 100):
    weights = 1./errs
    min_points = 2 
    linear_points = 15
    if len(mjds) > 1:
        max_diff = max(np.diff(mjds))
    else:
        max_diff = 0
    if len(mjds) > linear_points and max_diff < 20.:
        m_spline = splrep(mjds, mags, w = weights) #s smoothing factor?
        return m_spline
    elif len(mjds) > linear_points and max_diff >= 20.:
        m_spline = splrep(mjds, mags, w = weights, k=1)
        return m_spline
    elif len(mjds) <= linear_points and len(mjds) > min_points:
        m_spline = splrep(mjds, mags, w = weights, k=1)
        return m_spline
    else:
        return None


def generate_photometry_for_epoch_OSC(spec, valid_bands):
    phot = spec.light_curves
    phot_dict = {}
    for band in valid_bands:
        mjds, mags, errs = get_band_data(phot, band)
        m_spline = interp_LC(mjds, mags, errs)
        if mjds[0] < spec.mjd  and mjds[-1] > spec.mjd:
            m_smooth = splev(spec.mjd, m_spline)
            phot_dict[band] = float(m_smooth)
            #TODO: error estimation
        else:
            phot_dict[band] = np.nan
    return phot_dict

def generate_photometry_for_epoch(spec, valid_bands):
    phot = spec.homog_light_curves
    phot_dict = {}
    for band in valid_bands:
        if band in phot.keys():
            if band == 'V':
                band = 'V0'
            mjds, mags, errs = np.asarray(phot[band][0]), np.asarray(phot[band][1]), np.asarray(phot[band][2])
            if band == 'V0':
                band = 'V'
            m_spline = interp_LC(mjds, mags, errs)
            if mjds[0] < spec.mjd  and mjds[-1] > spec.mjd:
                m_smooth = splev(spec.mjd, m_spline)
                phot_dict[band] = float(m_smooth)
                #TODO: error estimation
            else:
                phot_dict[band] = np.nan
    return phot_dict


def get_filter_flux(wave, flux, filt):
    # sys.stdout = open(os.devnull, 'w')
    flux = filt.get_flux(wave, flux, axis = -1)
    # sys.stdout = sys.__stdout__
    return flux

def total_phot_offset(scale, spec, filters, target_mags):
    

    diff_sum = 0.
    for f in filters:
        if target_mags[f] is not None:
            filt_flux = get_filter_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], filters[f])
            # filt_flux = filters[f].get_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], axis = -1)
            mag_from_spec = -2.5 * np.log10(filt_flux) - filters[f].Vega_zero_mag # will this work for everything?
            diff = (mag_from_spec - target_mags[f])**2.
            diff_sum += diff

    

    return diff_sum


def valid_bands(spec):
    valid_bands = []
    band_waves = {"U": [3000.,4200.], 
                  "B": [3600.,5600.], 
                  "V": [4700.,7000.], 
                  "V0": [4700.,7000.],
                  "R": [5500.,9000.], 
                  "I": [7000.,9200.],
                  "u": [3000.,4000.], 
                  "g": [3600.,5600.], 
                  "r": [5300.,7000.], 
                  "i": [6600.,8400.], 
                  "z": [8000.,10500.]}
    min_wave = spec.wavelength[spec.x1]
    max_wave = spec.wavelength[spec.x2]

    for b in band_waves:
        if min_wave < band_waves[b][0] and max_wave > band_waves[b][1] and b in spec.light_curves.keys():
            valid_bands.append(b)


    return valid_bands


def scale_flux_to_photometry(spec, valid_bands):
    lib = pyphot.get_library()
    band_dict = {'U': 'GROUND_JOHNSON_U', 'B': 'GROUND_JOHNSON_B', 'V': 'GROUND_JOHNSON_V',
                 'R': 'GROUND_COUSINS_R', 'I': 'GROUND_COUSINS_I', 
                 'u': 'SDSS_u', 'g': 'SDSS_g', 'r': 'SDSS_r', 'i': 'SDSS_i', 'z': 'SDSS_z'}
    
    if len(valid_bands) > 0:
        mags_from_phot = generate_photometry_for_epoch(spec, valid_bands)
        
        valid_mjds = []
        for band in mags_from_phot:
            if ~np.isnan(mags_from_phot[band]):
                valid_mjds.append(band)
        if len(valid_mjds) > 0:
            filts = {}
            for b in valid_mjds:
                filts[b] = lib[band_dict.get(b)]

            guess = 1.e-14
            scale = opt.minimize(total_phot_offset, guess, args = (spec, filts, mags_from_phot), method = 'Nelder-Mead').x
            scale = scale[0]

            spec_phots = {}
            for b in valid_mjds:
                # final_flux = filts[b].get_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], axis = -1)
                final_flux = get_filter_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], filts[b])
                final = -2.5 * np.log10(final_flux) - filts[b].Vega_zero_mag
                spec_phots[b] = final
            print scale, spec_phots, mags_from_phot #for testing output
        else:
            scale = np.nan
            print scale, "mjd of spectrum outside photometric coverage"
    else:
        scale = np.nan
        mags_from_phot = None

    return scale, mags_from_phot

def MW_correction(spec, undo=False):
    #TODO:finish, apply before scaling

    wave = spec.wavelength[spec.x1:spec.x2]
    flux = spec.flux[spec.x1:spec.x2]
    ivar = spec.ivar[spec.x1:spec.x2]
    av_mw = spec.av_mw

    wave_u = wave*u.Angstrom
    flux_u = flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(wave_u, flux_u)
    spec1d_ivar = Spectrum1D.from_array(wave_u, ivar)
    red = ex.reddening(spec1d.wavelength, a_v=av_mw, r_v=3.1, model='f99')
    if not undo:
        flux_new = spec1d.flux*red
        ivar_new = spec1d_ivar.flux*(1./(red**2.))
    else:
        flux_new = spec1d.flux/red
        ivar_new = spec1d_ivar.flux/(1./(red**2.))

    spec.flux[spec.x1:spec.x2] = flux_new
    spec.ivar[spec.x1:spec.x2] = ivar_new

    return spec

def host_correction(spec, undo=False):
    #TODO:finish, apply before scaling

    wave = spec.wavelength[spec.x1:spec.x2]
    flux = spec.flux[spec.x1:spec.x2]
    ivar = spec.ivar[spec.x1:spec.x2]
    if spec.av_25 != None:
        Av_host = spec.av_25
        rv = 2.5
    elif spec.av_mlcs31 != None:
        Av_host = spec.av_mlcs31
        rv = 3.1
    elif spec.av_mlcs17 != None:
        Av_host = spec.av_mlcs17
        rv = 1.7
    else:
        print 'No Host Correction'
        return spec

    wave_u = wave*u.Angstrom
    flux_u = flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(wave_u, flux_u)
    spec1d_ivar = Spectrum1D.from_array(wave_u, ivar)
    red = ex.reddening(spec1d.wavelength, a_v=Av_host, r_v=rv, model='f99')
    if not undo:
        flux_new = spec1d.flux*red
        ivar_new = spec1d_ivar.flux*(1./(red**2.))
    else:
        flux_new = spec1d.flux/red
        ivar_new = spec1d_ivar.flux/(1./(red**2.))

    spec.flux[spec.x1:spec.x2] = flux_new
    spec.ivar[spec.x1:spec.x2] = ivar_new

    return spec


def make_colorbar(spec_array):
    params = []
    for spec in spec_array:
        # if spec.scale_to_phot != None:
        params.append(spec.phase)
    # print np.min(params), np.max(params)
    norm = matplotlib.colors.Normalize(vmin=np.min(params),vmax=np.max(params))
    c_m = matplotlib.cm.gist_rainbow
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    return s_m

def plot_light_curves_OSC(phot, fit = False, filt_list = None):

    kplot.basic_format()

    if filt_list is None:
        filt_list = phot.keys()

    for band in filt_list:
        mjds, mags, errs = get_band_data(phot, band)
        if fit and len(mjds) > 1:
            tnew = np.linspace(mjds[0], mjds[-1], 2000)
            m_spline = interp_LC(mjds, mags, errs)
            if m_spline is not None:
                m_smooth = splev(tnew, m_spline)
                plt.plot(tnew, m_smooth, color = color_dict[band])
        plt.errorbar(mjds, mags, yerr = errs, label = band, fmt='o', markersize=10, color = color_dict[band])
    plt.gca().invert_yaxis()
    plt.xlabel('MJD (days)', fontsize = 35)
    plt.ylabel('Magnitude', fontsize = 35)
    # plt.legend()
    plt.show()

def plot_light_curves(phot, name, spec_dates = None, fit = False, filt_list = None):
    
    kplot.basic_format()

    if filt_list is None:
        filt_list = phot.keys()

    for band in filt_list:
        mjds, mags, errs = np.asarray(phot[band][0]), np.asarray(phot[band][1]), np.asarray(phot[band][2])
        if fit and len(mjds) > 1:
            tnew = np.linspace(mjds[0], mjds[-1], 2000)
            m_spline = fc.interp_LC(mjds, mags, errs)
            if m_spline is not None:
                m_smooth = splev(tnew, m_spline)
                plt.plot(tnew, m_smooth, color = color_dict[band])
        plt.errorbar(mjds, mags, yerr = errs, label = band, fmt='o', markersize=10, color = color_dict[band])
    if spec_dates:
        for i, d in enumerate(spec_dates):
            if i == 0:
                plt.axvline(x=d, color = 'k', alpha = .7, label='Spectrum Epoch')
            else:
                plt.axvline(x=d, color = 'k', alpha = .7)
    plt.gca().invert_yaxis()
    plt.xlabel('MJD (days)', fontsize = 35)
    plt.ylabel('Magnitude', fontsize = 35)
    plt.legend(fontsize=20)
    plt.title(name, fontsize=20)
    plt.show()

def plot_calibrated_spectra(spec_array):

    kplot.basic_format()

    s_m = make_colorbar(spec_array)
    for i, spec in enumerate(spec_array):
        if spec.scale_to_phot:
            plt.plot(spec.wavelength[spec.x1:spec.x2], 
                     spec.scale_to_phot*spec.flux[spec.x1:spec.x2], 
                     color = s_m.to_rgba(spec.phase))
        # plt.plot(spec.wavelength[spec.x1:spec.x2], 
        #          spec.scale_to_phot*spec.flux[spec.x1:spec.x2])
    plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize = 35)
    plt.ylabel('Flux', fontsize = 35)
    plt.colorbar(s_m, label='Phase')
    plt.show()

def plot_spectra(spec_array):

    kplot.basic_format()

    s_m = make_colorbar(spec_array)
    for i, spec in enumerate(spec_array):
        plt.plot(spec.wavelength[spec.x1:spec.x2], spec.flux[spec.x1:spec.x2], 
                 color = s_m.to_rgba(spec.phase))
    plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize = 35)
    plt.ylabel('Flux', fontsize = 35)
    plt.colorbar(s_m, label='Phase')
    plt.show()


if __name__ == "__main__":

    query = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN where Spectra.SN = '2008r'"
    spec_array = kpora.grab(query, multi_epoch = True, make_corr = True, verbose=True, db_file = '../data/kaepora_v1.db')
    spec_array_dered = kpora.host_dereddening(spec_array, verbose=False, cutoff=2.)

    
    for spec in spec_array_dered:
        vbs = valid_bands(spec)
        scale, mags_from_phot = scale_flux_to_photometry(spec, vbs)
        spec.scale_to_phot = scale

    plot_calibrated_spectra(spec_array_dered)
