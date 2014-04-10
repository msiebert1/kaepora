import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import datafidelity as df

def tflag(wave, flux):

    tspec = pyfits.open('bstarr.fits')

    crval = tspec[0].header['CRVAL1']
    delta = tspec[0].header['CDELT1']
    tflux = tspec[0].data

    start = crval - 0.5*len(tflux)*delta
    stop = crval + 0.5*len(tflux)*delta
    twave = [(start+delta*i) for i in range(len(tflux))]
    twave = np.array(wave)

    spline_rep = interpolate.splrep(twave, tflux)
    new_telluric_flux = interpolate.splev(wave, spline_rep)

    #new_flux = df.gsmooth(wave_array, flux_array, var_y, vexp = 0.001)

    # Scaling telluric for comparison
    #telluric_flux = new_telluric_flux*np.median(flux_array)

SN1 = np.loadtxt('sn2007oo-20071107.34-fast.flm')

tflag(SN1[:,0], SN1[:, 1], SN1[:, 2])



    
