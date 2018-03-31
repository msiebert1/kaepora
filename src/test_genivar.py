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
    old_wave = old_wave*u.Angstrom        # wavelengths
    old_flux = old_flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(old_wave, old_flux)
    test_flux = test_dered.dered(sne, sn_name, spec1d.wavelength, spec1d.flux)  # Deredenning (see if sne in extinction files match the SN name)
#     new_flux = host_correction(sne, sn_name, old_wave, new_flux)

    # new_flux = old_flux
    new_flux = test_flux.value
    old_wave = old_wave.value
    old_wave = old_wave/(1.+z)

    old_flux = np.asarray(old_flux)
    new_flux = np.asarray(new_flux)
    s = scale_composites_in_range(old_flux, new_flux)
    old_flux = s*old_flux

    # plt.rc('font', family='serif')
    # fig, ax = plt.subplots(1,1)
    # fig.set_size_inches(10, 8, forward = True)
    # ax.get_yaxis().set_ticks([])
    # plt.plot(old_wave, old_flux, linewidth = 2, color = 'r')
    # plt.plot(old_wave, new_flux, linewidth = 2, color = '#3F5D7D')
    # plt.ylabel('Relative Flux')
    # plt.xlabel('Wavelength ' + "($\mathrm{\AA}$)")
    # # plt.savefig('../../Paper_Drafts/MW_corr.png', dpi = 300, bbox_inches = 'tight')
    # plt.show()

    av = .1294 #2006sr
    # av = 2.9711 #2005a
    name = '2006sr'
    # name = '2005a'
    host_wave = old_wave*u.Angstrom        # wavelengths
    host_flux = new_flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(host_wave, host_flux)
    new_flux_host, new_ivar_host = test_dered.host_correction(av, 2.5, name, spec1d.wavelength, spec1d.flux, [0])

    new_flux = np.asarray(new_flux)
    new_flux_host = np.asarray(new_flux_host.value)
    s = scale_composites_in_range(new_flux_host, new_flux)
    new_flux_host = s*new_flux_host

    norm = 1./np.amax(new_flux_host)
    new_flux_host = new_flux_host*norm
    new_flux = new_flux*norm
    old_flux = old_flux*norm

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
    plt.plot(old_wave, old_flux, linewidth = 2, color = '#d95f02')
    plt.plot(old_wave, new_flux, linewidth = 2, color = '#1b9e77')
    plt.plot(host_wave.value, new_flux_host, linewidth = 2, color = '#7570b3')
    plt.ylabel('Relative Flux', fontsize = 30)
    plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
    plt.savefig('../../Paper_Drafts/red_corr.pdf', dpi = 300, bbox_inches = 'tight')
    # plt.ylim([-.2,1.01])
    # plt.savefig('../../Paper_Drafts/red_corr_large.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()

    # new_wave = old_wave/(1.+z)  # Deredshifting
    new_wave = old_wave
    new_error = old_error  # Placeholder if it needs to be changed
    norm = 1./np.amax(new_flux)
    new_flux = new_flux*norm
    new_ivar = df.genivar(new_wave, new_flux, new_error)  # generate new inverse variance
    #var = new_flux*0+1
    newdata = prep.Interpo(new_wave, new_flux, new_ivar)  # Do the interpolation
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
    plt.plot(old_wave, new_flux, linewidth = 2, color = 'r')
    plt.plot(newdata[0], newdata[1], linewidth = 2, color = '#3F5D7D')
    plt.ylabel('Relative Flux', fontsize = 30)
    plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
    plt.savefig('../../Paper_Drafts/interp.pdf', dpi = 300, bbox_inches = 'tight')
    # plt.ylim([-.3,1.])
    # plt.savefig('../../Paper_Drafts/interp_large.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()

#    print 'new spectra',newdata
    return newdata, snr

#reddening / interpolation plots
c = 299792.458
#sn2006sr-20061220.097-ui.flm
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

# fname = '../data/spectra/bsnip/sn2005a-20050115.357-br.flm'
# spectrum = np.loadtxt(fname)
# # sn_name = '2006sr'
# sn_name = '2005a'
# source = 'bsnip'
# bsnip_vals = read_bsnip_data('obj_info_table.txt')
# if is_number(fname.split('-')[1][:8]):
#     longdate = fname.split('-')[1][:8]
# else:
#     longdate = fname.split('-')[2][:8]
# data = bsnip_vals[sn_name.lower()+'-'+longdate]
# redshift = data[1]/c
# newdata, snr = compprep(spectrum, sn_name, redshift, source)

raise TypeError

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
plt.plot(old_wave, real_var, linewidth = 2, color = '#7570b3', alpha = .8)
plt.plot(old_wave, new_var, linewidth = 6, color = 'k')
plt.ylabel('Variance', fontsize = 30)
plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
# plt.savefig('../../Paper_Drafts/genvar.pdf', dpi = 300, bbox_inches = 'tight')
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