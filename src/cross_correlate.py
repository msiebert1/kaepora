import composite as comp 
import matplotlib.pyplot as plt
import numpy as np
import pyfits
from math import *
from scipy import interpolate

#not working correctly
def sky_correlation():
	test_query = " "
	query = "SELECT * FROM Supernovae " + test_query

	SN_Array = comp.grab(query)

	sky = pyfits.open('../personal/AdamSnyder/kecksky.fits')
	crval = sky[0].header['CRVAL1']
	delta = sky[0].header['CDELT1']
	skyflux = sky[0].data[0]
	start = crval
	stop = crval + ceil(len(skyflux)*delta)
	skywave = [(start+delta*i) for i in range(len(skyflux))]
	# Find wavelength overlap
	spline_rep = interpolate.splrep(skywave, skyflux)

	telluric_spec = np.loadtxt('../data/etc/telluric_spec.dat')
	telluric_spec = np.transpose(telluric_spec)
	tell_wave = telluric_spec[0]
	tell_flux = telluric_spec[1]
	spline_rep_tell = interpolate.splrep(skywave, skyflux)

	corr_dict = {}
	for SN in  SN_Array:
		good = np.where((SN.wavelength >= 7550.) & (SN.wavelength <= 7650.))
		# good = np.where((SN.wavelength >= skywave[0]) & (SN.wavelength <= skywave[-1]))

		sky_flux_interp = interpolate.splev(SN.wavelength[good], spline_rep_tell)

		# SN.sky_corr = np.corrcoef(SN.flux[good]*1./np.average(SN.flux[good]), sky_flux_interp)[0][1]
		SN.tell_corr = np.corrcoef(SN.flux[good]*1./np.average(SN.flux[good]), sky_flux_interp)[0][1]
		print SN.tell_corr
	for SN in SN_Array:
		if np.absolute(SN.tell_corr) > .4:
			plt.plot(SN.wavelength, SN.flux)
			plt.show()
sky_correlation()