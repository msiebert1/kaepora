import datafidelity as df
import matplotlib.pyplot as plt
import numpy as np 

def find_extrema(SN):
	wavelength = SN.wavelength[SN.x1:SN.x2]
	flux = SN.flux[SN.x1:SN.x2]
	varflux = 1./SN.ivar[SN.x1:SN.x2]
	sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
	# sm_flux = flux

	flux_deriv = np.gradient(sm_flux, wavelength)
	flux_deriv_abs = np.absolute(flux_deriv)
	zero_inds = np.where(flux_deriv_abs < 1e-8)[0]

	plt.plot(wavelength[zero_inds], sm_flux[zero_inds],'go', markersize=25.)
	plt.plot(wavelength,flux,'r-', linewidth=4.)
	plt.plot(wavelength,sm_flux, linewidth=4.)
	plt.show()

	plt.plot(wavelength,flux_deriv, linewidth=4.)
	plt.show()

def measure_si_velocity(SN):
	wavelength = SN.wavelength[SN.x1:SN.x2]
	flux = SN.flux[SN.x1:SN.x2]
	varflux = 1./SN.ivar[SN.x1:SN.x2]
	sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
	si_range = np.where((wavelength > 5850.) & (wavelength < 6400.))
	si_wave = wavelength[si_range]
	si_flux = sm_flux[si_range]
	si_min = np.amin(si_flux)
	si_min_index = np.where(si_flux == si_min)
	si_min_wave = si_wave[si_min_index][0]

	c = 299792. # km/s
	si_rest_wave = 6355. #Angstroms

	v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
	# plt.plot(wavelength, flux)
	# plt.plot(wavelength, sm_flux)
	# plt.plot(si_min_wave, si_min, 'o', color='orange')
	# plt.show()
	return v

def measure_weak_si_velocity(SN):
	wavelength = SN.wavelength[SN.x1:SN.x2]
	flux = SN.flux[SN.x1:SN.x2]
	varflux = 1./SN.ivar[SN.x1:SN.x2]
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
	# plt.plot(wavelength, flux)
	# plt.plot(wavelength, sm_flux)
	# plt.plot(si_min_wave, si_min, 'o')
	# plt.show()
	return v

	#carbon II is 6580A (2006bt good example)

def measure_si_velocity_from_raw(wavelength,flux,z):
	# varflux = 1./SN.ivar[SN.x1:SN.x2]
	# plt.plot(wavelength,flux)
	wavelength = wavelength/(1.+z)
	varflux = np.ones(len(flux), float)
	sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
	si_range = np.where((wavelength > 5830.) & (wavelength < 6230.))
	si_wave = wavelength[si_range]
	si_flux = sm_flux[si_range]
	si_min = np.amin(si_flux)
	si_min_index = np.where(si_flux == si_min)
	si_min_wave = si_wave[si_min_index][0]

	c = 299792. # km/s
	si_rest_wave = 6355. #Angstroms
	v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
	# v = c*(si_min_wave - si_rest_wave)/si_rest_wave
	plt.plot(wavelength, flux)
	plt.plot(wavelength, sm_flux)
	plt.plot(si_min_wave, si_min, 'o', color='orange')
	plt.show()
	return v, si_min_wave

def measure_C_velocity_from_raw(wavelength,flux,z):
	# varflux = 1./SN.ivar[SN.x1:SN.x2]
	# plt.plot(wavelength,flux)
	wavelength = wavelength/(1.+z)
	varflux = np.ones(len(flux), float)
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
	plt.plot(wavelength, flux)
	plt.plot(wavelength, sm_flux)
	plt.plot(C_min_wave, C_min, 'o', color='orange')
	plt.show()
	return v, C_min_wave

spec_file = r'..\..\Observing\Keck\01_15_18\2018ha_comb.dat'
with open(spec_file) as spec:
	spectrum = np.transpose(np.loadtxt(spec))
	print spectrum
	wavelength = spectrum[0]
	flux = spectrum[1]
	plt.plot(wavelength, flux)
	plt.show()
	# z = 0.005274
	# # print  measure_si_velocity_from_raw(wavelength, flux, z)
	# print  measure_C_velocity_from_raw(wavelength, flux, z)