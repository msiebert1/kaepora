import datafidelity as df 
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import os
import prep as prep
import pandas as pd
from astropy import units as u
from specutils import extinction as ex
from specutils import Spectrum1D
import test_dered
import math
import scipy.interpolate as inter

def scale_composites_in_range(data, comp):
    scales = []
    guess = 1.
    s = opt.minimize(sq_residuals_in_range, guess, args = (data, comp), 
                 method = 'Nelder-Mead').x
    return s

def sq_residuals_in_range(s, data, comp):
    data = s*data
    res = data - comp
    sq_res = res*res
    return np.sum(sq_res)

def read_bsnip_data(data_file):
    """
    Returns a dict keyed on 'snname-longdate' ie. 1989a-19890427
    """
    data1 = pd.read_fwf('../data/info_files/obj_info_table.txt', names=('SN name', 'Type',
                                        'Host Galaxy', 'Host Morphology', 'Redshift', 'Reddening',
                                        'Discovery Date'), colspecs=((0, 10), (10, 17), (17, 51),
                                        (51, 58), (58, 63), (63, 69), (69, 97)))

    data2 = pd.read_fwf('../data/info_files/spec_info_table.txt',
                        names=('SN name', 'Year', 'Month', 'Day',
                                       'MJD of Spectrum', 'Phase of Spectrum',),
                        colspecs=((0, 9), (10, 15), (15, 18), (18, 25),
                                         (25, 36), (36, 43)))
    dmat1 = data1.as_matrix()
    dmat2 = data2.as_matrix()
    dmat2 = dmat2.tolist()
    #format month as '04' instead of '4'
    #do the same for day (5.000 -> 05.000)
    for line in dmat2:
        if line[2] < 10:
            line[2] = ''.join(['0', str(line[2])])
        if line[3] < 10:
            line[3] = ''.join(['0', str(line[3])])
    for line in dmat2:
        fulldate = ''.join([str(line[1]), str(line[2]), str(line[3])[:2]])
        line.append(fulldate)

    bsnip_dict_rs = {}
    bsnip_dict = {}
    #split matrix lines to populate dicts
    for line in dmat1:
        key = line[0].split()[1].lower()
        rs = line[4]
        bsnip_dict_rs[key] = rs
    for line in dmat2:
        key1 = line[0].split()[1].lower()
        key2 = line[len(line)-1]
        full_key = key1+'-'+key2
        # bsnip_dict[full_key] = [key1, bsnip_dict_rs[key1], line[len(line)-2]]
        bsnip_dict[full_key] = [key1, bsnip_dict_rs[key1], line[len(line)-2], line[len(line)-3]]
    return bsnip_dict

#test more spectra
# scales = []
# i = 0

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def clip(wave, flux, ivar):
    # Create an array of all ones
    var = np.ones(len(flux), float)
    
    # Create 2 smoothed fluxes, of varying vexp
    sflux = df.gsmooth(wave, flux, var, 0.002)
    
    # Take the difference of the two fluxes and smooth
    err = abs(flux - sflux)
    serr = df.gsmooth(wave, err, var, 0.008)

    # Find the wavelengths that need to be clipped (omitting 5800-6000 region)
    index = np.where(((err/serr > 5.5) & (wave < 5800.0)) | ((err/serr > 5.5) & (wave > 6000.0)))
    index_na = np.where((err/serr > 10.5) & ((wave > 5800.0) & (wave < 6000.0)))
    bad_wave = wave[index]
    bad_wave_na = wave[index_na]

    # Find indices for general clipping
    # bad = np.array([], int)
    bad_ranges = [] # if don't need it to be a numpy array (A.S.)
    
    for i in range(len(bad_wave)):
        # bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))
        bad_ranges.append((bad_wave[i]-8, bad_wave[i]+8)) # instead save each point as tuple of desired range (A.S.)

    for i in range(len(bad_wave_na)):
        # bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))
        bad_ranges.append((bad_wave_na[i]-8, bad_wave_na[i]+8))

    # Set ivar to 0 for those points and return
#    ivar[bad] = 0
#    return ivar

    plt.plot(wave, flux)
    for wave_tuple in bad_ranges:
#        print wave_tuple
        clip_points = np.where((wave > wave_tuple[0]) & (wave < wave_tuple[1]))
        #make sure not at edge of spectrum
        flux[clip_points] = np.interp(wave[clip_points], [wave_tuple[0], wave_tuple[1]], [flux[clip_points[0][0]-1], flux[clip_points[0][-1]+1]])
        #deweight data (but not to 0), somewhat arbitrary
        ivar[clip_points] = np.interp(wave[clip_points], [wave_tuple[0], wave_tuple[1]], [ivar[clip_points[0][0]-1], ivar[clip_points[0][-1]+1]])

    plt.plot(wave, flux)
    plt.show()
    return wave, flux, ivar # return bad_ranges instead of setting ivar[bad] = 0 (A.S.)

def Interpo (wave, flux, ivar, sn_name, plot=False):
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

    #ivar = clip(wave, flux, ivar) #clip bad points in flux (if before interpolation)
    # ivar = clipmore(wave,flux,ivar)    
    # bad_points = clip(wave, flux, ivar)  # if returned bad points range instead of ivar
#    print 'ivar', ivar
#    print 'bad points', bad_points
    #ivar[ivar < 0] = 0 # make sure no negative points

    good_data = np.where((wave >= lower) & (wave <= upper))  #creates an array of wavelength values between minimum and maximum wavelengths from new spectrum

    influx = inter.splrep(wave[good_data], flux[good_data])  # creates b-spline from new spectrum

    inivar = inter.splrep(wave[good_data], ivar[good_data])  # doing the same with the inverse varinces

    # extrapolating returns edge values
    inter_flux = inter.splev(wavelength, influx, ext = 3)    # fits b-spline over wavelength range
    inter_ivar = inter.splev(wavelength, inivar, ext = 3)   # doing the same with errors

    scale = 10.
    if plot:
        missing_data = np.where((wavelength < lower) | (wavelength > upper))
        inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
        inter_ivar[missing_data] = float('NaN')
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
        ax.xaxis.set_ticks(np.arange(np.round(wave[0],-3),np.round(wave[-1],-3),1000))
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
        plt.plot(wave, scale*flux, linewidth = 2, color = 'r', label='Before Interpolation')
        plt.plot(wavelength, scale*inter_flux, linewidth = 2, color = '#3F5D7D', label='After Interpolation')
        plt.ylabel('Relative Flux', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([wave[0]-200,wave[-1]+200])
        if sn_name == '2005a':
            plt.ylim([-.3*scale,1.05*scale])
            plt.legend(loc=4, fontsize=20)
            # plt.savefig('../../../Paper_Drafts/reprocessing/interp_large_av.pdf', dpi = 300, bbox_inches = 'tight')
        else:
            plt.ylim([-.05*scale,1.05*scale])
            plt.legend(loc=1, fontsize=20)
            # plt.savefig('../../../Paper_Drafts/reprocessing/interp_small_av.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()

#    inter_ivar = clip(wavelength, inter_flux, inter_var) #clip bad points (if do after interpolation)

    # Then the below code (or something similar) would do it (A.S.)
    inter_flux_old = np.copy(inter_flux)

#     for wave_tuple in bad_points:
# #        print wave_tuple
#         zero_points = np.where((wavelength > wave_tuple[0]) & (wavelength < wave_tuple[1]))
#         #make sure not at edge of spectrum
#         inter_flux[zero_points] = np.interp(wavelength[zero_points], [wave_tuple[0], wave_tuple[1]], [inter_flux[zero_points[0][0]-1], inter_flux[zero_points[0][-1]+1]])
#         #deweight data (but not to 0), somewhat arbitrary
#         inter_ivar[zero_points] = np.interp(wavelength[zero_points], [wave_tuple[0], wave_tuple[1]], [inter_ivar[zero_points[0][0]-1], inter_ivar[zero_points[0][-1]+1]])
#         # inter_ivar[zero_points] = 0

    inter_ivar[inter_ivar < 0] = 0  # make sure there are no negative points!
    
#    place = np.where((wavelength > 5800.0 ) & (wavelength < 6000.0 ))
#    print inter_ivar[place]
    missing_data = np.where((wavelength < lower) | (wavelength > upper))
    inter_flux[missing_data] = float('NaN')  # set the bad values to NaN !!!
    inter_ivar[missing_data] = float('NaN')

    if plot:
        plt.rc('font', family='serif')
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches(10, 8, forward = True)
        plt.minorticks_on()
        plt.xticks(fontsize = 20)
        ax.xaxis.set_ticks(np.arange(np.round(wavelength[0],-3),np.round(wavelength[-1],-3),1000))
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
        plt.plot(wavelength, scale*inter_flux_old, linewidth = 2, color = 'r', label='Before Clipping')
        plt.plot(wavelength, scale*inter_flux, linewidth = 2, color = '#3F5D7D', label='After Clipping')
        plt.ylabel('Relative Flux', fontsize = 30)
        plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
        plt.xlim([lower-200, upper+200])
        if sn_name == '2005a':
            plt.ylim([-.3*scale,1.05*scale])
            plt.legend(loc=4, fontsize=20)
            # plt.savefig('../../../Paper_Drafts/reprocessing/clip_large_av.pdf', dpi = 300, bbox_inches = 'tight')
        else:
            plt.ylim([-.05*scale,1.05*scale])
            plt.legend(loc=1, fontsize=20)
            # plt.savefig('../../../Paper_Drafts/reprocessing/clip_small_av.pdf', dpi = 300, bbox_inches = 'tight')
        plt.show()

#    print inter_ivar[place]

    output = np.array([wavelength, inter_flux, inter_ivar])  # put the interpolated data into the new table

    return output  # return new table

def compprep(spectrum, sn_name, z, source):
    old_wave = spectrum[:, 0]       # wavelengths
    old_flux = spectrum[:, 1]   # fluxes
    try:
        old_error = spectrum[:, 2]  # check if supernovae has error array
    except IndexError:
        old_error = np.array([0])  # if not, set default
    old_ivar = df.genivar(old_wave, old_flux, old_error)  # generate inverse variance
    snr = prep.getsnr(old_flux, old_ivar)

    if source == 'cfa':  # choosing source dataset
#        z = ReadParam()
        sne = prep.ReadExtin('extinction.dat')
    if source == 'bsnip':
        sne = prep.ReadExtin('extinctionbsnip.dat')
    if source == 'csp':
        sne = prep.ReadExtin('extinctioncsp.dat')
        old_wave *= 1+float(z)  # Redshift back
    if source == 'uv':
        sne = prep.ReadExtin('extinctionuv.dat')
    if source == 'other':
        sne = prep.ReadExtin('extinctionother.dat')

#     host_reddened = ReadExtin('../data/info_files/ryan_av.txt')
    newdata = []
    dered_wave = old_wave/(1.+z)
    new_error = old_error  # Placeholder if it needs to be changed
    norm = 1./np.amax(old_flux)
    new_flux = old_flux*norm
    new_ivar = df.genivar(dered_wave, new_flux, new_error)  # generate new inverse variance
    #var = new_flux*0+1
    dered_wave, new_flux, new_ivar = clip(dered_wave, new_flux, new_ivar)
    newdata = Interpo(dered_wave, new_flux, new_ivar, sn_name, plot=True)  # Do the interpolation

    interp_wave = newdata[0]*u.Angstrom        # wavelengths
    interp_flux = newdata[1]*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(interp_wave, interp_flux)
    test_flux = test_dered.dered(sne, sn_name, spec1d.wavelength, spec1d.flux)  # Deredenning (see if sne in extinction files match the SN name)
#     new_flux = host_correction(sne, sn_name, old_wave, new_flux)

    # new_flux = old_flux
    mw_flux = test_flux.value
    interp_wave = interp_wave.value
    interp_flux = interp_flux.value

    interp_flux = np.nan_to_num(interp_flux)
    mw_flux = np.nan_to_num(mw_flux)
    s = scale_composites_in_range(mw_flux, interp_flux)
    mw_flux = s*mw_flux

    if sn_name == '2006sr':
        av = .1294 #2006sr
        name = '2006sr'
    else:
        av = 2.9711 #2005a
        name = '2005a'
    # name = '2005a'
    host_wave = interp_wave*u.Angstrom        # wavelengths
    host_flux = mw_flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(host_wave, host_flux)
    new_flux_host, new_ivar_host = test_dered.host_correction(av, 2.5, name, spec1d.wavelength, spec1d.flux, [0])

    new_flux_host = np.asarray(new_flux_host.value)
    s = scale_composites_in_range(new_flux_host, interp_flux)
    new_flux_host = s*new_flux_host

    scale = 10.
    valid_data = np.where(interp_wave > 3600)
    norm = 1./np.nanmax(new_flux_host[valid_data])
    new_flux_host = new_flux_host*norm*scale
    mw_flux = mw_flux*norm*scale
    interp_flux = interp_flux*norm*scale

    interp_flux[interp_flux == 0] = np.nan
    mw_flux[mw_flux == 0] = np.nan
    new_flux_host[new_flux_host == 0] = np.nan

    plt.rc('font', family='serif')
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(10, 8, forward = True)
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    ax.xaxis.set_ticks(np.arange(np.round(interp_wave[0],-3),np.round(interp_wave[-1],-3),1000))
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
    plt.plot(interp_wave, interp_flux, linewidth = 2, color = '#d95f02', label = 'Before Dereddening')
    plt.plot(interp_wave, mw_flux, linewidth = 2, color = '#1b9e77', label = 'MW Corrected')
    plt.plot(interp_wave, new_flux_host, linewidth = 2, color = '#7570b3', label = 'Host Corrected')
    plt.ylabel('Relative Flux', fontsize = 30)
    plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
    # plt.savefig('../../Paper_Drafts/red_corr.pdf', dpi = 300, bbox_inches = 'tight')
    plt.xlim([dered_wave[0]-200,dered_wave[-1]+200])
    plt.legend(loc=1, fontsize=20)
    if sn_name == '2005a':
        plt.ylim([-.05*scale,1.05*scale])
        # plt.savefig('../../../Paper_Drafts/reprocessing/red_corr_large_av.pdf', dpi = 300, bbox_inches = 'tight')
    else:
        plt.ylim([-.05*scale,1.05*scale])
        # plt.savefig('../../../Paper_Drafts/reprocessing/red_corr_small_av.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()

    # new_wave = old_wave/(1.+z)  # Deredshifting

#    print 'new spectra',newdata
    return newdata, snr

#reddening / interpolation plots
c = 299792.458
#sn2006sr-20061220.097-ui.flm


###CODE FOR INTERP AND REDDENING PLOTS
fname = '../data/spectra/bsnip/sn2006sr-20061220.097-ui.flm'
spectrum = np.loadtxt(fname)
sn_name = '2006sr'
source = 'bsnip'
bsnip_vals = read_bsnip_data('obj_info_table.txt')
if is_number(fname.split('-')[1][:8]):
    longdate = fname.split('-')[1][:8]
else:
    longdate = fname.split('-')[2][:8]
data = bsnip_vals[sn_name.lower()+'-'+longdate]
redshift = data[1]/c
newdata, snr = compprep(spectrum, sn_name, redshift, source)

fname = '../data/spectra/bsnip/sn2005a-20050115.357-br.flm'
spectrum = np.loadtxt(fname)
# sn_name = '2006sr'
sn_name = '2005a'
source = 'bsnip'
bsnip_vals = read_bsnip_data('obj_info_table.txt')
if is_number(fname.split('-')[1][:8]):
    longdate = fname.split('-')[1][:8]
else:
    longdate = fname.split('-')[2][:8]
data = bsnip_vals[sn_name.lower()+'-'+longdate]
redshift = data[1]/c
newdata, snr = compprep(spectrum, sn_name, redshift, source)
raise TypeError
### ENDS HERE



fname = '../data/spectra/cfa/sn2003kg/sn2003kg-20031129.12-fast.flm'
spectrum = np.loadtxt(fname)

old_wave = spectrum[:, 0]
old_flux = spectrum[:, 1]

real_error = spectrum[:, 2]
if real_error[0] != 0.:
    old_error = np.zeros(len(old_wave), float)
    new_ivar = df.genivar(old_wave, old_flux, old_error)
    new_var = 1./new_ivar
    real_var = real_error*real_error
    scale = scale_composites_in_range(new_var, real_var)
    # new_var = new_var*scale
    new_var = new_var*2.02
    print scale

norm = 1./np.amax(real_var)
real_var = real_var*norm
new_var = new_var*norm
plt.rc('font', family='serif')
fig, ax = plt.subplots(1,1)
fig.set_size_inches(10, 8, forward = True)
plt.minorticks_on()
plt.xticks(fontsize = 20)
ax.xaxis.set_ticks(np.arange(np.round(old_wave[0],-3),np.round(old_wave[-1],-3),1000))
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
plt.plot(old_wave, real_var, linewidth = 2, color = '#7570b3', alpha = .8, label = 'CfA')
plt.plot(old_wave, new_var, linewidth = 6, color = 'k', label='Generated')
plt.ylabel('Variance', fontsize = 30)
plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
plt.xlim([old_wave[0]-200.,old_wave[-1]+200.])
plt.ylim([0.,1.1])
plt.legend(fontsize=20)
# plt.savefig('../../../Paper_Drafts/reprocessing/genvar.pdf', dpi = 300, bbox_inches = 'tight')
plt.show()

# plt.subplot(2,1,1)
# plt.plot(output[0], output[1])
# plt.subplot(2,1,2)
# plt.plot(output[0], output[2])
# plt.plot()
# plt.show()

# for path, subdirs, files in os.walk('../data/spectra/cfa'):
# 	if i < 1000:
# 		for fname in files:

# 			fname = os.path.join(path, fname)
# 			if fname.endswith('.flm'):
# 				spectrum = np.loadtxt(fname)

# 				old_wave = spectrum[:, 0]
# 				old_flux = spectrum[:, 1]

# 				real_error = spectrum[:, 2]
# 				if real_error[0] != 0.:
# 					old_error = np.zeros(len(old_wave), float)
# 					new_ivar = df.genivar(old_wave, old_flux, old_error)

# 					new_var = 1./new_ivar
# 					real_var = real_error*real_error
# 					scale = scale_composites_in_range(new_var, real_var)
# 					new_var = new_var*scale
# 					print scale
# 					scales.append(scale)
# 					i += 1
# 	else:
# 		break

			# plt.subplot(2,1,1)
			# plt.plot(old_wave, old_flux)
			# plt.subplot(2,1,2)
			# plt.plot(old_wave, real_var)
			# plt.plot(old_wave, new_var)
			# plt.show()


# print np.average(scales)
# print np.std(scales)

#avg: 2.02 std: 1.055