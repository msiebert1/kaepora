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

def compprep(spectrum, sn_name, z, source):
    old_wave = spectrum[:, 0]       # wavelengths
    old_flux = spectrum[:, 1]   # fluxes
    try:
        old_error = spectrum[:, 2]  # check if supernovae has error array
    except IndexError:
        old_error = np.array([0])  # if not, set default
    if sn_name == '2011fe':
        old_error = np.sqrt(old_error)
    old_ivar = df.genivar(old_wave, old_flux, old_error)  # generate inverse variance
    snr = prep.getsnr(old_flux, old_ivar)

    if source == 'cfa':  # choosing source dataset
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

    newdata = []
    old_wave = old_wave*u.Angstrom
    old_flux = old_flux*u.Unit('W m-2 angstrom-1 sr-1')
    spec1d = Spectrum1D.from_array(old_wave, old_flux)
    test_flux = test_dered.dered(sne, sn_name, spec1d.wavelength, spec1d.flux)
    new_flux = test_flux.value
    old_wave = old_wave.value
    old_wave = old_wave/(1.+z)

    old_flux = np.asarray(old_flux)
    new_flux = np.asarray(new_flux)
    # s = scale_composites_in_range(old_flux, new_flux)
    # old_flux = s*old_flux
    new_wave = old_wave/(1.+z)
    new_error = old_error
    new_ivar = df.genivar(new_wave, new_flux, new_error)
    newdata = prep.Interpo(new_wave, new_flux, new_ivar)
    return newdata, snr

def build_redshift_dict(bsnipdict, cfadict):
    """
    Creates a dictionary of the form snname: redshift from the cfa and bsnip data
    """
    rsd = {}
    for item in cfadict:
        rsd[item] = float(cfadict[item][0])
    for item in bsnipdict:
        if item not in rsd:
            rsd[item] = float(bsnipdict[item])
            rsd[item.lower()] = float(bsnipdict[item])
    return rsd

def read_cfa_info(data_file, dates_file):
    """
    Reads Cfa SN information from separate files.
    Output dict format:
    # key,     0         1           2       3        4        5       6     7        8        9
    # SN,  zhel,  tmax(B),  +/-  ref.,  Dm15, +/-  ref.,  M_B   +/-,   B-V,
     10        11       12       13              14
    +/-,  Bm-Vm, +/-,   Phot. ref   Obs Date
    """
    with open(dates_file) as dates:
        lines = dates.readlines()
        date_dict = {}
        for line in lines:
            if not line.startswith('#'):
                cleandate = line.split()

                if cleandate:
                    # date_dict[cleandate[0]] = cleandate[1]
                    date_dict[cleandate[0]] = cleandate[1]

    with open(data_file) as data:
        lines = data.readlines()
        sndict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                sndict[sndata[0].lower()] = sndata[1:]

    return sndict, date_dict


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

def create_short_bsnip_dict(bsnipdict):
    newdict = {}
    for k, v in bsnipdict.iteritems():
        newdict[k.split('-')[0]] = v[1]/c
    return newdict

def find_SN(fname, source=None, csplist=None):
    """
    Returns SN name, either from cspdict if source is a csp spectra
    or from slicing the file name if source is Cfa or bsnip
    """
    if source == 'csp':
        snname = csplist[0]
        return snname[2:]
    elif source == 'other':
        snname = fname.replace('_', '-').split('-')
        print snname
        if snname[0][:2] == 'sn':
            return snname[0][2:].lower()
        else:
            return snname[0].lower()
    else:
        snname = fname.replace('_', '-').split('-')
        if snname[0][:3] == 'snf':
            namelist = [snname[0], snname[1]]
            snname = '-'.join(namelist).upper()
        else:
            snname = snname[0][2:]

        # return snname
        return snname.lower()

def build_NED_redshift_dict(NED_file):
    with open(NED_file) as data:
        lines = data.readlines()
        NED_red_dict = {}
        for line in lines:
            sn_name = line.split()[0]
            redshift = float(line.split()[1])
            NED_red_dict[sn_name] = redshift

    return NED_red_dict

c = 299792.458
sndict, date_dict = read_cfa_info('../data/spectra/cfa/cfasnIa_param.dat',
                                  '../data/spectra/cfa/cfasnIa_mjdspec.dat')
bsnip_vals = read_bsnip_data('obj_info_table.txt')
short_bsnip_dict = create_short_bsnip_dict(bsnip_vals)
rsd = build_redshift_dict(short_bsnip_dict, sndict)
NED_red_dict = build_NED_redshift_dict('../data/info_files/NED_redshift_info.txt')
#sn2006sr-20061220.097-ui.flm
fname = '../data/spectra/other/sn2011fe-20110929-snifs.dat'
spectra = np.loadtxt(fname)
sn_name = find_SN('sn2011fe-20110929-snifs.dat')
redshift = None
if sn_name in rsd:
    redshift = rsd[sn_name]
print sn_name, redshift
if sn_name in NED_red_dict:
    redshift = NED_red_dict[sn_name]
print sn_name, redshift
source = 'uv'
newdata, snr = compprep(spectra, sn_name, redshift, source)
print snr
wavelength = newdata[0,:]
flux = newdata[1,:]
ivar = newdata[2,:]
plt.plot(wavelength,flux)
plt.show()
plt.plot(wavelength,ivar)
plt.show()