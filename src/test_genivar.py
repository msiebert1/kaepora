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
import pyfits
from scipy import interpolate
from math import *

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
# c = 299792.458
# #sn2006sr-20061220.097-ui.flm
# fname = '../data/spectra/bsnip/sn2006sr-20061220.097-ui.flm'
# spectrum = np.loadtxt(fname)
# sn_name = '2006sr'
# source = 'bsnip'
# bsnip_vals = read_bsnip_data('obj_info_table.txt')
# if is_number(fname.split('-')[1][:8]):
#     longdate = fname.split('-')[1][:8]
# else:
#     longdate = fname.split('-')[2][:8]
# data = bsnip_vals[sn_name.lower()+'-'+longdate]
# redshift = data[1]/c
# newdata, snr = compprep(spectrum, sn_name, redshift, source)

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

# raise TypeError

def gsmooth(x_array, y_array, var_y, vexp , nsig = 5.0):
    
    # Check for zero variance points, and set to 1E-20
    for i in range(len(var_y)):
        if var_y[i] == 0:
            var_y[i] = 1E-20
            # var_y[i] = 1E-31
    
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
    scale = 50.*med_error #fudge factor provides reasonable scaling of sky lines
    add_flux = scale*add_flux

    # Add sky flux to the error
    new_error = error
    new_error[good] = error[good] + add_flux

    return new_error

def genivar(wavelength, flux, varflux , vexp = 0.0008, nsig = 5.0): 

    # Check to see if it has a variance already
    if len(np.where(varflux == 0)[0]) > 0:
        varflux = np.ones(len(wavelength))
    else:
        ivar = 1 / (varflux**2)
        return ivar
    
    # Smooth original flux
    # new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)
    new_flux = gsmooth(wavelength, flux, varflux, .002, nsig)
    
    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)
    
    # Smooth noise to find the variance
    # sm_error = gsmooth(wavelength, error, varflux, vexp, nsig)
    sm_error = gsmooth(wavelength, error, varflux, .007, nsig)

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
    
    sm_var_new = 2.02*sm_error*sm_error #multiply by ideal fudge factor to account for consistently low variance estimate

    # Inverse variance
    ivar = 1./(sm_var_new)
    
    # Return generated variance
    return ivar

def getsnr(flux, ivar):
    sqvar = map(math.sqrt, ivar)
    snr = flux/(np.divide(1.0, sqvar))
    snr_med = np.median(snr)
    return snr_med

fname = '../data/spectra/bsnip/sn2001fg-20011118.251-joined.flm'
spectrum = np.loadtxt(fname)
bad_files = [u'sn2001eu-20011022.342-joined.flm', u'sn2001fg-20011118.251-joined.flm', u'sn2002hd-20021109-joined.flm', u'sn2002he-20021106.649-joined.flm', u'sn2008ha-20081205-ui.flm', 
             u'sn1994D-19940602.18-fast.flm', u'sn1995ak-19951222.23-fast.flm', u'sn1995bd-19951225.27-fast.flm', u'sn1996bo-19961107.15-fast.flm', u'sn1996by-19970108.33-fast.flm', 
             u'sn2000fa-20001225.42-fast.flm', u'sn2001C-20010325.19-fast.flm', u'sn2001C-20010328.18-fast.flm', u'sn2001E-20010124.48-fast.flm', u'sn2001E-20010225.37-fast.flm', 
             u'sn2002dl-20020618.44-fast.flm', u'sn2002do-20020710.33-fast.flm', u'sn2002eu-20020905.49-fast.flm', u'sn2002fb-20020929.42-fast.flm', u'sn2002fk-20021112.36-fast.flm', 
             u'sn2002hw-20021114.10-fast.flm', u'sn2002kf-20030113.31-fast.flm', u'sn2003ag-20030209.41-fast.flm', u'sn2003ch-20030428.15-fast.flm', u'sn2003gj-20030704.40-mmt.flm', 
             u'sn2003hu-20030928.12-fast.flm', u'sn2003it-20031123.14-fast.flm', u'sn2003kf-20040216.18-fast.flm', u'sn2003kf-20040318.12-fast.flm', u'sn2003kz-20031219.51-fast.flm', 
             u'sn2003S-20030127.49-fast.flm', u'sn2003Y-20030131.35-fast.flm', u'sn2004bp-20040511.21-fast.flm', u'sn2004dt-20040908.49-fast.flm', u'sn2004dt-20041015.34-fast.flm', 
             u'sn2004dt-20041212.18-fast.flm', u'sn2004ef-20040915.30-fast.flm', u'sn2004ef-20040922.22-fast.flm', u'sn2004gs-20050106.40-fast.flm', u'sn2004H-20040127.46-fast.flm', 
             u'sn2005am-20050517.14-fast.flm', u'sn2005cc-20050617.18-fast.flm', u'sn2005cc-20050628.25-fast.flm', u'sn2005cc-20050709.21-fast.flm', u'sn2005M-20050209.33-fast.flm', 
             u'sn2005M-20050409.19-fast.flm', u'sn2005mz-20060123.15-fast.flm', u'sn2005na-20060130.24-fast.flm', u'sn2006ac-20060226.34-fast.flm', u'sn2006az-20060403.30-fast.flm', 
             u'sn2006bb-20060402.20-fast.flm', u'sn2006bq-20060427.35-fast.flm', u'sn2006br-20060430.18-fast.flm', u'sn2006cf-20060521.27-fast.flm', u'sn2006cf-20060522.26-fast.flm', 
             u'sn2006cj-20060523.33-fast.flm', u'sn2006cj-20060601.28-fast.flm', u'sn2006dt-20060721.20-fast.flm', u'sn2006dv-20060725.43-fast.flm', u'sn2006H-20060203.12-fast.flm', 
             u'sn2006lf-20061221.21-fast.flm', u'sn2006lf-20061222.24-fast.flm', u'sn2006R-20060129.53-fast.flm', u'sn2006te-20070109.49-fast.flm', u'sn2007af-20070423.32-fast.flm', 
             u'sn2007al-20070314.26-fast.flm', u'sn2007bj-20070609.28-fast.flm', u'sn2007bj-20070613.29-fast.flm', u'sn2007ci-20070614.18-fast.flm', u'sn2007F-20070116.54-fast.flm', 
             u'sn2007F-20070121.55-fast.flm', u'sn2007F-20070309.39-fast.flm', u'sn2007F-20070415.31-fast.flm', u'sn2007if-20071010.31-fast.flm', u'sn2007kk-20071015.32-fast.flm', 
             u'sn2008A-20080226.14-fast.flm', u'sn2008ae-20080210.36-fast.flm', u'sn2008ae-20080214.32-fast.flm', u'sn2008at-20080305.30-fast.flm', u'sn2008C-20080402.14-mmt.flm', 
             u'sn2008E-20080203.40-fast.flm', u'sn2008Y-20080209.41-fast.flm', u'sn2008Y-20080229.35-fast.flm', u'sn2008Z-20080416.16-fast.flm', u'2005cf_20050601_3243_9720_00.dat', 
             u'sn2004dt-20040820-hst.flm', u'sn2004dt-20040823-hst.flm', u'sn2004ef-20040914-hst.flm', u'sn2004ef-20040918-hst.flm', u'sn2005cf-20050603-hst.flm', 
             u'sn2005cf-20050605-hst.flm', u'sn2005cf-20050607-hst.flm', u'sn2005cf-20050611-hst.flm', u'sn2005cf-20050614-hst.flm', u'sn2005m-20050128-hst.flm', u'sn2005m-20050131-hst.flm']
old_wave = spectrum[:, 0]
old_flux = spectrum[:, 1]

real_error = spectrum[:, 2]
# if real_error[0] != 0.:
old_error = np.zeros(len(old_wave), float)
new_ivar = genivar(old_wave, old_flux, old_error)
new_var = 1./new_ivar
real_var = real_error*real_error
# scale = scale_composites_in_range(new_var, real_var)
# new_var = new_var*scale
new_var = new_var*2.02
# print scale

# norm = 1./np.nanmax(real_var)
# real_var = real_var*norm
# new_var = new_var*norm
plt.plot(old_wave, real_var)
plt.plot(old_wave, new_var)
plt.show()
# plt.plot(old_wave, new_var, linewidth = 6, color = 'k')
# plt.plot(old_wave, real_var, linewidth = 2, color = '#7570b3', alpha = .8)
# # plt.savefig('../../Paper_Drafts/genvar.pdf', dpi = 300, bbox_inches = 'tight')
# plt.show()

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