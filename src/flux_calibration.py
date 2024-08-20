import matplotlib
import matplotlib.pyplot as plt
import kaepora as kpora
import kaepora_plot as kplot
import composite
import flux_calibration as fc
# import pyphot
import pysynphot
import numpy as np
from scipy.interpolate import splrep, splev
from scipy.interpolate import interp2d
import scipy
import scipy.optimize as opt
# from specutils import extinction as ex
from specutils import Spectrum1D
from dust_extinction.parameter_averages import F99
from astropy import units as u
import sys,os
import copy


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
              'Y': 'pink',

              'W1': 'gray',
              'W2': 'black',
              'M2': 'darkgray',

              'C': 'limegreen',
              'Js': 'gold',
              'Ks': 'lightgray'}

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

# color_dict = {b'U': 'magenta',
#           b'B': 'blue',
#           b'V': 'green',
#           b'R': 'red',
#           b'I': 'brown',

#           b'u': 'purple',
#           b'g': 'seagreen',
#           b'r': 'crimson',
#           b'i': 'orangered',
#           b'z': 'darkred',

#           b"u'": 'purple',
#           b"g'": 'seagreen',
#           b"r'": 'crimson',
#           b"i'": 'orangered',
#           b"z'": 'darkred',

#           b'H': 'yellow',
#           b'K': 'teal',
#           b'J': 'crimson',
#           b'V': 'orange',
#           b'V0': 'orange',
#           b'Y': 'pink',

#           b'W1': 'gray',
#           b'W2': 'black',
#           b'M2': 'darkgray',

#           b'C': 'limegreen',
#           b'Js': 'gold',
#           b'Ks': 'lightgray',
#           b'Jrc2': 'lavender',
#           b'Hdw': 'darkblue',
#           b'Ydw': 'darkgreen',
#           b'Jdw': 'darkgray',
#              }

def get_band_data(phot, band):
    mjds = []
    mags = []
    errs = []

    telescope_dict = {"None": 0}
    system_dict = {"None": 0}
    observatory_dict = {"None": 0}
    instrument_dict = {"None": 0}

    for mag in phot[band][1]:
        if mag[1].get(b'telescope', None) is not None:
            tscope = mag[1].get(b'telescope', None)
            if tscope not in telescope_dict:
                telescope_dict[tscope] = 1
            else:
                telescope_dict[tscope] = telescope_dict[tscope] + 1
        else:
            telescope_dict['None'] = telescope_dict['None'] + 1

        if mag[1].get(b'system', None) is not None:
            sys = mag[1].get(b'system', None)
            if sys not in system_dict:
                system_dict[sys] = 1
            else:
                system_dict[sys] = system_dict[sys] + 1
        else:
            system_dict['None'] = system_dict['None'] + 1

        if mag[1].get(b'observatory', None) is not None:
            obs = mag[1].get(b'observatory', None)
            if obs not in observatory_dict:
                observatory_dict[obs] = 1
            else:
                observatory_dict[obs] = observatory_dict[obs] + 1
        else:
            observatory_dict['None'] = observatory_dict['None'] + 1

        if mag[1].get(b'instrument', None) is not None:
            inst = mag[1].get(b'instrument', None)
            if inst not in instrument_dict:
                instrument_dict[inst] = 1
            else:
                instrument_dict[inst] = instrument_dict[inst] + 1
        else:
            instrument_dict['None'] = instrument_dict['None'] + 1
    print (band)
    print (telescope_dict)
    print (system_dict)
    print (observatory_dict)
    print (instrument_dict)

    max_sys = max(system_dict, key=lambda k: system_dict[k])
    print (max_sys)

    for i, mag in enumerate(phot[band][1]):
        sys = mag[1].get('system', "None")
        if  sys == max_sys:
            mjds.append(float(phot[band][0][i]))
            mags.append(float(mag[0]))
            if b'e_magnitude' in mag[1]:
                errs.append(float(mag[1][b'e_magnitude']))
            elif b'upperlimit' in mag[1] and mag[1][b'upperlimit']:
                errs.append(np.nan)
            else:
                errs.append(0.5)
    print (len(mjds))

    mjds = np.asarray(mjds)
    mags = np.asarray(mags)
    errs = np.asarray(errs)

    detection_inds = np.where(~np.isnan(errs))[0]
    mjds = mjds[detection_inds]
    mags = mags[detection_inds]
    errs = errs[detection_inds]
    print
    return mjds, mags, errs


def interp_LC(mjds, mags, errs, s = 100, linear=False):
    weights = 1./errs
    if not linear:
        m_spline = splrep(mjds, mags, w = weights)
    else:
        m_spline = splrep(mjds, mags, w = weights, k=1)
    return m_spline


def generate_photometry_for_epoch_OSC(spec, valid_bands, template=None):
    phot = spec.light_curves
    phot_dict = {}
    for band in valid_bands:
        mjds, mags, errs = get_band_data(phot, band)
        if template == 'salt' and spec.mjd > 48:
            m_spline_ext = interp_LC(mjds[-4:], mags[-4:], errs[-4:],linear=True)
        else:
            m_spline = interp_LC(mjds, mags, errs)

        if mjds[0] < spec.mjd  and mjds[-1] > spec.mjd:
            m_smooth = splev(spec.mjd, m_spline)
            phot_dict[band] = float(m_smooth)
            #TODO: error estimation
        elif mjds[-1] < spec.mjd and template =='salt':
            m_ext = splev(spec.mjd, m_spline_ext, ext=0)
            phot_dict[band] = float(m_ext)
        else:
            phot_dict[band] = np.nan
    return phot_dict

def generate_photometry_for_epoch(spec, valid_bands, template=None):
    # phot = spec.homog_light_curves
    phot = spec.event_data["Homogenized_Photometry"]
    phot_dict = {}
    for band in valid_bands:
        if band in phot.keys():
            if band == 'V':
                band = 'V0'
            mjds, mags, errs = np.asarray(phot[band][0]), np.asarray(phot[band][1]), np.asarray(phot[band][2])
            if band == 'V0':
                band = 'V'

            if template == 'salt' and spec.mjd > 49:
                m_spline_ext = interp_LC(mjds[-4:], mags[-4:], errs[-4:],linear=True)
            else:
                m_spline = interp_LC(mjds, mags, errs)

            # print (mjds[0], mjds[-1], spec.mjd)
            if mjds[0] < spec.mjd  and mjds[-1] > spec.mjd:
                m_smooth = splev(spec.mjd, m_spline)
                phot_dict[band] = float(m_smooth)
                #TODO: error estimation
            elif mjds[-1] < spec.mjd and template =='salt':
                m_ext = splev(spec.mjd, m_spline_ext)
                phot_dict[band] = float(m_ext)
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

def total_phot_offset_pysyn(scale, spec, filters, target_mags):
    
    diff_sum = 0.

    if scale < 0:
        scale = 1.e-14
    sp = pysynphot.ArraySpectrum(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], fluxunits='flam')

    for f in filters:
        bp = filters[f]
        name = bp.name.split('/')[-1].split('_')[0] + ','+ bp.name.split('/')[-1].split('_')[1]
        if target_mags[f] is not None:
            spec_obs = pysynphot.Observation(sp, bp, force='extrap')
            if 'johnson' in name:
                mag_from_spec = spec_obs.effstim('vegamag')
            elif 'sdss' in name:
                mag_from_spec = spec_obs.effstim('abmag')

            diff = (mag_from_spec - target_mags[f])**2.
            diff_sum += diff

    return diff_sum


def valid_bands(spec):
    valid_bands = []
    # band_waves = {"U": [3000.,4200.], 
    #               "B": [3600.,5600.], 
    #               "V": [4700.,7000.], 
    #               "V0": [4700.,7000.],
    #               "R": [5500.,9000.], 
    #               "I": [7000.,9200.],
    #               "u": [3000.,4000.], 
    #               "g": [3600.,5600.], 
    #               "r": [5300.,7000.], 
    #               "i": [6600.,8400.], 
    #               "z": [8000.,10500.]}
    band_waves = {b"U": [3000.,4200.], 
                  b"B": [3600.,5600.], 
                  b"V": [4700.,7000.], 
                  b"V0": [4700.,7000.],
                  b"R": [5500.,9000.], 
                  b"I": [7000.,9200.],
                  b"u": [3000.,4000.], 
                  b"g": [3600.,5600.], 
                  b"r": [5300.,7000.], 
                  b"i": [6600.,8400.], 
                  b"z": [8000.,10500.]}
    min_wave = spec.wavelength[spec.x1]
    max_wave = spec.wavelength[spec.x2]

    for b in band_waves:
        if min_wave < band_waves[b][0] and max_wave > band_waves[b][1] and b in spec.event_data['Homogenized_Photometry'].keys():
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
            print (scale, spec_phots, mags_from_phot) #for testing output
        else:
            scale = np.nan
            print (scale, "mjd of spectrum outside photometric coverage")
    else:
        scale = np.nan
        mags_from_phot = None

    return scale, mags_from_phot


def scale_flux_to_photometry_pysyn(spec, valid_bands, template=None):
    # lib = pyphot.get_library()
    # band_dict = {'U': 'johnson,u', 'B': 'johnson,b', 'V': 'johnson,v',  'V0': 'johnson,v',
    #              'R': 'johnson,r', 'I': 'johnson,i', 
    #              'u': 'sdss,u', 'g': 'sdss,g', 'r': 'sdss,r', 'i': 'sdss,i', 'z': 'sdss,z'}
    band_dict = {b'U': 'johnson,u', b'B': 'johnson,b', b'V': 'johnson,v',  b'V0': 'johnson,v',
                 b'R': 'johnson,r', b'I': 'johnson,i', 
                 b'u': 'sdss,u', b'g': 'sdss,g', b'r': 'sdss,r', b'i': 'sdss,i', b'z': 'sdss,z'}
    
    if len(valid_bands) > 0:
        mags_from_phot = generate_photometry_for_epoch(spec, valid_bands, template=template)
        
        valid_mjds = []
        for band in mags_from_phot:
            if ~np.isnan(mags_from_phot[band]):
                valid_mjds.append(band)
        if len(valid_mjds) > 0:
            filts = {}
            for b in valid_mjds:
                # filts[b] = lib[band_dict.get(b)]
                filts[b] =  pysynphot.ObsBandpass(band_dict.get(b))

            guess = 1.e-14
            scale = opt.minimize(total_phot_offset_pysyn, guess, args = (spec, filts, mags_from_phot), method = 'Nelder-Mead').x
            scale = scale[0]

            spec_phots = {}
            # for f in filts:
            #     bp = filters[f]
            #     name = bp.name.split('/')[-1].split('_')[0] + ','+ bp.name.split('/')[-1].split('_')[1]
            #     if target_mags[f] is not None:
            #         spec_obs = pysynphot.Observation(sp, bp, force='extrap')
            #         if 'johnson' in name:
            #             mag_from_spec = spec_obs.effstim('vegamag')
            #         elif 'sdss' in name:
            #             mag_from_spec = spec_obs.effstim('abmag')

            #     final_flux = filts[b].get_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], axis = -1)
            #     final_flux = get_filter_flux(spec.wavelength[spec.x1:spec.x2], scale*spec.flux[spec.x1:spec.x2], filts[b])
            #     final = -2.5 * np.log10(final_flux) - filts[b].Vega_zero_mag
            #     spec_phots[b] = final
            # print scale, spec_phots, mags_from_phot #for testing output
        else:
            scale = np.nan
            print (scale, "mjd of spectrum outside photometric coverage")
    else:
        scale = np.nan
        mags_from_phot = None

    return scale, mags_from_phot


def MW_correction(spec, undo=False):
    #TODO:finish, apply before scaling

    wave = spec.wavelength[spec.x1:spec.x2]
    flux = spec.flux[spec.x1:spec.x2]
    ivar = spec.ivar[spec.x1:spec.x2]
    av_mw = spec.event_data['Av_MW']

    wave_u = wave*u.Angstrom
    flux_u = flux*u.Unit('W m-2 angstrom-1 sr-1')
    # spec1d = Spectrum1D.from_array(wave_u, flux_u)
    # spec1d_ivar = Spectrum1D(spectral_axis=wave, flux=flux)
    # red = ex.reddening(spec1d.wavelength, a_v=av_mw, r_v=3.1, model='f99')

    ext = F99(Rv=3.1)
    red = ext.extinguish(wave*u.AA, Av = av_mw)

    if not undo:
        # flux_new = spec1d.flux*red
        # ivar_new = spec1d_ivar.flux*(1./(red**2.))
        flux_new = flux*red
        ivar_new = ivar*(1./(red**2.))
    else:
        # flux_new = spec1d.flux/red
        # ivar_new = spec1d_ivar.flux/(1./(red**2.))
        flux_new = flux/red
        ivar_new = ivar/(1./(red**2.))

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
        print ('No Host Correction')
        return spec

    wave_u = wave*u.Angstrom
    flux_u = flux*u.Unit('W m-2 angstrom-1 sr-1')
    # red = ex.reddening(spec1d.wavelength, a_v=Av_host, r_v=rv, model='f99')

    ext = F99(Rv=3.1)
    red = ext.extinguish(wave*u.AA, Av = Av_host)

    if not undo:
        # flux_new = spec1d.flux*red
        # ivar_new = spec1d_ivar.flux*(1./(red**2.))
        flux_new = flux*red
        ivar_new = ivar*(1./(red**2.))
    else:
        # flux_new = spec1d.flux/red
        # ivar_new = spec1d_ivar.flux/(1./(red**2.))
        flux_new = flux/red
        ivar_new = ivar/(1./(red**2.))

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

def fit_lc_lecho(data, fit = False, filt_list = None, linear=False):

    kplot.basic_format()

    if filt_list is None:
        filt_list = phot.keys()

    splines = []
    for i, band in enumerate(filt_list):
        mjds, mags, errs = data[i]
        if fit and len(mjds) > 1:
            tnew = np.linspace(mjds[0], mjds[-1], 2000)
            m_spline = interp_LC(mjds, mags, errs, linear=linear)
            if m_spline is not None:
                m_smooth = splev(tnew, m_spline)
                plt.plot(tnew, m_smooth, color = color_dict[band])
                splines.append(m_spline)
        plt.errorbar(mjds, mags, yerr = errs, label = band, fmt='o', markersize=10, color = color_dict[band])
    plt.gca().invert_yaxis()
    plt.xlabel('MJD (days)', fontsize = 35)
    plt.ylabel('Magnitude', fontsize = 35)
    # plt.xlim(55800,55840)
    # plt.ylim(14,9.9)
    plt.show()
    return splines


def plot_light_curves(phot, name, spec_dates = None, fit = False, filt_list = None):
    
    kplot.basic_format()

    if filt_list is None:
        filt_list = phot.keys()

    for band in filt_list:
        mjds, mags, errs = np.asarray(phot[band][0]), np.asarray(phot[band][1]), np.asarray(phot[band][2])
        if fit and len(mjds) > 1:
            tnew = np.linspace(mjds[0], mjds[-1], 2000)
            tnew_ext = np.linspace(mjds[-3], 100, 2000)
            m_spline = fc.interp_LC(mjds, mags, errs)
            m_spline_ext = fc.interp_LC(mjds[-4:], mags[-4:], errs[-4:],linear=True)
            if m_spline is not None:
                m_smooth = splev(tnew, m_spline)
                m_ext = splev(tnew_ext, m_spline_ext,ext=0)
                plt.plot(tnew, m_smooth, color = color_dict[band])
                plt.plot(tnew_ext, m_ext, color = 'k')
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


def create_salt_phot(x1):
    phase_s, wavelength_salt, flux_salt = np.loadtxt('../data/info_files/SALT/salt2-4/salt2_template_0.dat', unpack = True)
    phase_s_corr, wavelength_salt_corr, salt_corr = np.loadtxt('../data/info_files/SALT/salt2-4/salt2_template_1.dat', unpack = True)
    flux_salt = flux_salt + x1*salt_corr #dm15=1.1

    salt_temps = {}
    phase_list = []
    first = True
    for i, phase in enumerate(phase_s):
        if phase not in phase_list:
            phase_list.append(phase)
            if not first:
                salt_temps[float(phase_s[i-1])] = [temp_wave, temp_flux]
            temp_wave = [wavelength_salt[i]]
            temp_flux = [flux_salt[i]]
        else:
            temp_wave.append(wavelength_salt[i])
            temp_flux.append(flux_salt[i])
        first = False
        
    i=-19.
    phases_s = []

    # band_dict = {
    #     b'B':pysynphot.ObsBandpass('johnson,b'),
    #     b'V':pysynphot.ObsBandpass('johnson,v'),
    #     b'R':pysynphot.ObsBandpass('johnson,r'),
    #     b'I':pysynphot.ObsBandpass('johnson,i') # this is wacky right now
    # }

    band_dict = {
        b'B':pysynphot.ObsBandpass('johnson,b'),
        b'V':pysynphot.ObsBandpass('johnson,v'),
        b'R':pysynphot.ObsBandpass('johnson,r')
    }

    salt_lcs = {}
    while i <= np.amax(list(salt_temps.keys())):
        sp = pysynphot.ArraySpectrum(np.asarray(salt_temps[i][0]), np.asarray(np.asarray(salt_temps[i][1])), fluxunits='flam')
        for band in band_dict.keys():
            bp = band_dict[band]
            spec_obs = pysynphot.Observation(sp, bp, force='extrap')
            mag_from_spec = spec_obs.effstim('vegamag')
            if salt_lcs.get(band, None) is not None:
                salt_lcs[band].append(mag_from_spec)
            else:
                salt_lcs[band] = [mag_from_spec]

        phases_s.append(i)
        i+=1
        
    phot = {}
    err_h = 0.01*np.ones(len(phases_s))
    for band in band_dict.keys():
        phot[band] = [phases_s,salt_lcs[band],err_h,'salt']
    return phot


def create_hsiao_phot():
    t_h, U_h, B_h, V_h, R_h, I_h, Y_h, J_h, H_h, K_h = np.loadtxt('../data/info_files/hsiao_lc_template.dat', unpack = True)
    err_h = 0.01*np.ones(len(t_h))
    phot = {}
    phot[b'B'] = [t_h[5:],B_h[5:],err_h[5:],'hsiao']
    phot[b'V'] = [t_h[5:],V_h[5:],err_h[5:],'hsiao']
    phot[b'R'] = [t_h[5:],R_h[5:],err_h[5:],'hsiao']
    phot[b'I'] = [t_h[5:],I_h[5:],err_h[5:],'hsiao']
    return phot

def scale_composites_to_photometry(composite_dict, verbose=False, template=None):
    spec_array = composite.prelim_norm(composite_dict.values())
    for spec in spec_array:
        spec.mjd = np.average(spec.phase_array[spec.x1:spec.x2])
        vbs = fc.valid_bands(spec)
        scale, mags_from_phot = fc.scale_flux_to_photometry_pysyn(spec, vbs, template=template)
        if ~np.isnan(scale):
            spec.scale_to_phot = scale
        else:
            spec.scale_to_phot = None

    for spec in spec_array:
        if verbose:
            print (spec.scale_to_phot)
        spec.flux *= spec.scale_to_phot
        spec.ivar /= spec.scale_to_phot**2

    return spec_array

def prepare_spec(spec_array,waverange=[3500,9000],plot=False):
    specnew = {}
    phases = []
    for spec in spec_array:
        snew = {}
        wave = spec.wavelength[spec.x1:spec.x2]
        flux = spec.flux[spec.x1:spec.x2]
        ivar = spec.ivar[spec.x1:spec.x2]
        if plot:
            plt.plot(wave,flux)
#         output, scale, var_final = prep.Interpo_flux_conserving(wave,flux,ivar,dw=dw)
#         wave,flux,ivar = output[0],output[1],output[2]
        
        idx = np.array((wave>waverange[0]) & (wave<waverange[1]))
        idx_nan = ~np.isnan(flux)
        
        idx_nan_w1 = (wave<spec.wavelength[spec.x1])
        flux[idx_nan_w1] = np.median(flux[idx_nan][0:10])
        print (np.median(flux[idx_nan][0:10]))
        ivar[idx_nan_w1] = ivar[idx_nan][0]/1e3
        
        idx_nan_w2 = (wave>spec.wavelength[spec.x2])
        flux[idx_nan_w2] = np.median(flux[idx_nan][-10:])
        print (np.median(flux[idx_nan][-10:]))
        ivar[idx_nan_w2] = ivar[idx_nan][-1]/1e3
        
        print (wave[idx])
        snew['wavelength_interp'] = wave[idx]
        snew['flux_interp'] = flux[idx]
        snew['ivar_interp'] = ivar[idx]
        
#         flux[idx_nan] = 0.
#         ivar[idx_nan] = 0.
        
        phase = np.average(spec.phase_array[spec.x1:spec.x2])
        phases.append(phase)
        
        specnew[phase] = snew
        if plot:
            plt.plot(wave[idx],flux[idx])
#             plt.plot(wave,flux)
        plt.title(str(phase))
        plt.show()
    return specnew, np.sort(phases), wave[idx]

def linear_interpolation(spec, waves, phases, pint_range=[-20,90]):
    nwave = np.shape(spec[list(spec.keys())[0]]['wavelength_interp'])[0]
    nphase = len(list(spec.keys()))
    arr2d = []
    for i,phase in enumerate(np.sort(list(spec.keys()))):
        # nan_locs = np.where(np.isnan(spec[phase]['flux_interp']))[0]
        # for nanl in nan_locs:
        #     if nanl - 5 < 0:
        #         spec[phase]['flux_interp'][nanl] = np.nanmedian(spec[phase]['flux_interp'][nanl:nanl+5])
        #     elif nanl + 5 >= len(spec[phase]['flux_interp']):
        #         spec[phase]['flux_interp'][nanl] = np.nanmedian(spec[phase]['flux_interp'][nanl-5:nanl])
        #     else: 
        #         spec[phase]['flux_interp'][nanl] = np.nanmedian(spec[phase]['flux_interp'][nanl-5:nanl+5])
        arr2d.append(spec[phase]['flux_interp'])

    weight2d = []
    for i,phase in enumerate(np.sort(list(spec.keys()))):
        weight2d.append(spec[phase]['weight_interp'])
    x = waves
    y = phases
    X, Y = np.meshgrid(x, y)
    Z = np.asarray(arr2d)
    W = np.asarray(weight2d)

    plt.figure(figsize=[10,10])
    plt.imshow(Z,aspect=20, origin='lower', interpolation='none')
    plt.show()
    # Z = np.sin(np.pi*X/2) * np.exp(Y/2)
    print (np.shape(x))
    print(np.shape(y))
    print (np.shape(Z))
    print (np.shape(W))

    x2 = waves
    y2 = np.arange(pint_range[0], pint_range[1], 1)
    # y2 = np.arange(0, 90, 1)

    # f = interp2d(x, y, Z, kind='linear')
    # Z2 = f(x2, y2)

    tx = x
    ty = y
    splrep = scipy.interpolate.LSQBivariateSpline(X.ravel(), Y.ravel(), Z.ravel(), tx, ty, w=W.ravel(), kx=1, ky=1)
    Z2 = np.transpose(splrep(x2,y2))

    X2, Y2 = np.meshgrid(x2, y2)

    return X2,Y2,Z2,X,Y,Z

def plot_spectra(spec_array):

    fig, ax, = kplot.basic_format()

    s_m = make_colorbar(spec_array)
    for i, spec in enumerate(spec_array):
        plt.plot(spec.wavelength[spec.x1:spec.x2], spec.flux[spec.x1:spec.x2], 
                 color = s_m.to_rgba(spec.phase))
    plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize = 35)
    plt.ylabel('Flux', fontsize = 35)
    cax=ax.inset_axes([1, 0., 0.05, 1])
    plt.colorbar(s_m, cax=cax,label='Phase')
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
