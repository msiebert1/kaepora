#!/usr/bin/env python
from astropy.io import fits
from optparse import OptionParser
import numpy as np
from stissplice import splicer, tools
import matplotlib.pyplot as plt
from astropy import units as u
import astropy.constants as const
from scipy import signal
import spectral_analysis as sa
import datafidelity as df

if __name__ == "__main__":

    description = "For converting fits to dat"
    usage = "%prog    \t [option] \n Recommended syntax: %prog"
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option("-z", "--redshift", dest="redshift", action="store_true",
                      help='redshift of target')
    option, args = parser.parse_args()
    _redshift= option.redshift

    fits_file = fits.open(args[0])
    x1d_head = fits_file[0].header
    data = fits_file['sci',1].data
    #JWST stuff
    # fits_file = fits.open(args[0])
    # x1d_head = fits_file[0].header
    # x1d_head = fits_file[1].header
    # # data = fits_file['sci',1].data
    # data = fits_file[1].data
    # wave = fits_file[9].data*1e10


    # data_og = data*u.watt/(u.meter**3)
    # data = data_og.to(u.erg/u.s/u.cm**2/u.Angstrom)
    # # data = data_og.to(u.jansky)
    # data_fnu = (wave**2/const.c)*data_og
    # # data = (wave**2/const.c)*data_og

    # if _redshift:
    #     wave = wave/(1.+float(args[1]))

    # # plt.figure(figsize=[15,8])

    # wave = wave[~np.isnan(data)]
    # data = data[~np.isnan(data)]

    # plt.figure(figsize=[15,8])
    # plt.plot(wave,data)
    # plt.show()

    # plt.plot(wave,data,drawstyle='steps-mid')
    # sdata = df.gsmooth(wave, data.value, None, 0.05)
    # plt.plot(wave,sdata,drawstyle='steps-mid')
    # # plt.plot(wave,data_fnu,drawstyle='steps-mid')
    # plt.savefig(args[0][0:-4]+'.png')
    # plt.show()

    # resids = np.absolute(data.value - sdata)
    # sig = np.std(resids[(wave>8000) & (wave<9000)])

    # tol=5
    # new_wave = wave[resids/sig < tol]
    # new_data = data[resids/sig < tol]
    # plt.figure(figsize=[15,8])
    # plt.plot(wave,data,drawstyle='steps-mid')
    # plt.plot(new_wave,new_data,drawstyle='steps-mid')
    # plt.show()
    # data = data_fnu
    # wave = new_wave
    # data = new_data

    if x1d_head.get('OPT_ELEM', None) is not None and x1d_head['OPT_ELEM'].startswith('E'):

        spectrum = splicer.read_spectrum(args[0].split('_x1d')[0],'./')
        spliced_spectrum = splicer.splice_pipeline(args[0].split('_x1d')[0],'./')

        # for s in spectrum:
        #     plt.plot(s['wavelength'], s['flux'], alpha=0.3)
        # _ = plt.xlabel(r'Wavelength (${\rm \AA}$)')
        # _ = plt.ylabel(r'Flux density (erg s$^{-1}$ cm$^{-2}$ ${\rm \AA}^{-1}$)')
        # plt.plot(spliced_spectrum['WAVELENGTH'], spliced_spectrum['FLUX'], color='k')
        # plt.show()

        wave = np.asarray(spliced_spectrum['WAVELENGTH'])
        flux = np.asarray(spliced_spectrum['FLUX'])
        spec = np.asarray([wave,flux])
        np.savetxt(args[0][0:-4] + 'dat', np.transpose(spec))

    else:
        if 'prism_clear' in args[0]:
            wave = wave[~np.isnan(data)]
            flux = data[~np.isnan(data)]
            wave = np.asarray(wave)
            flux = np.asarray(flux)
            spec = np.asarray([wave,flux])
        elif x1d_head.get('TELESCOP').strip() == 'IUE':
            wave = np.asarray(data[0][0])
            flux = np.asarray(data[0][1])
            spec = np.asarray([wave,flux])
        else:
            wave = np.asarray(data['WAVELENGTH'])[0]
            flux = np.asarray(data['FLUX'])[0]
            spec = np.asarray([wave,flux])
        np.savetxt(args[0][0:-4] + 'dat', np.transpose(spec))


