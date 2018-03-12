import datafidelity as df
import matplotlib.pyplot as plt
import numpy as np 
import copy

def find_extrema(wavelength,sm_flux):
	#input a smoothed spectrum
	dw = 2.
	deriv_flux = np.diff(sm_flux)/dw
	abs_deriv = np.absolute(deriv_flux)
	extrema_locs = np.where(abs_deriv < .0001)
	maxima = []
	for loc in extrema_locs[0]:
		if loc != 0 and loc != len(deriv_flux)-1:
			if deriv_flux[loc-1] > 0 and deriv_flux[loc+1] < 0:
				maxima.append(loc)

	return maxima


def measure_si_ratio(wavelength,flux, varflux = None, vexp=.002, smooth=True):
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)

	if smooth:
		sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
	else:
		sm_flux = flux

	maxima = find_extrema(wavelength,sm_flux)

	si_range_1 = np.where((wavelength > 5500.) & (wavelength < 5700.))
	si_range_2 = np.where((wavelength > 5700.) & (wavelength < 6000.))
	si_range_3 = np.where((wavelength > 6000.) & (wavelength < 6500.))

	# si_wave_1 = wavelength[si_range_1]
	# si_wave_2 = wavelength[si_range_2]
	# si_wave_3 = wavelength[si_range_3]

	# si_flux_1 = sm_flux[si_range_1]
	# si_flux_2 = sm_flux[si_range_2]
	# si_flux_3 = sm_flux[si_range_3]

	# si_max_1 = np.amax(si_flux_1)
	# si_max_2 = np.amax(si_flux_2)
	# si_max_3 = np.amax(si_flux_3)

	# si_max_index_1 = np.where(si_flux_1 == si_max_1)
	# si_max_index_2 = np.where(si_flux_2 == si_max_2)	
	# si_max_index_3 = np.where(si_flux_3 == si_max_3)

	# si_max_wave_1 = si_wave_1[si_max_index_1][0]
	# si_max_wave_2 = si_wave_2[si_max_index_2][0]
	# si_max_wave_3 = si_wave_3[si_max_index_3][0]

	m1 = None
	m2 = None
	m3 = None
	for m in maxima:
		if m in si_range_1[0]:
			m1 = m
		if m in si_range_2[0]:
			m2 = m
		if m in si_range_3[0]:
			m3 = m

	if m1 is None or m2 is None or m3 is None:
		sm_flux = df.gsmooth(wavelength, flux, varflux, vexp= vexp+.001)
		maxima = find_extrema(wavelength,sm_flux)
		for m in maxima:
			if m in si_range_1[0]:
				m1 = m
			if m in si_range_2[0]:
				m2 = m
			if m in si_range_3[0]:
				m3 = m
	if m1 is None or m2 is None or m3 is None:
		print "Could not find maximum in a specified range!"
		return np.nan

	si_max_wave_1 = wavelength[m1]
	si_max_wave_2 = wavelength[m2]
	si_max_wave_3 = wavelength[m3]

	si_max_1 = sm_flux[m1]
	si_max_2 = sm_flux[m2]
	si_max_3 = sm_flux[m3]

	weak_si_trough = np.where((wavelength >= si_max_wave_1) & (wavelength < si_max_wave_2))[0]
	strong_si_trough = np.where((wavelength >= si_max_wave_2) & (wavelength < si_max_wave_3))[0]

	interp_flux = copy.copy(flux)
	interp_flux[weak_si_trough] = np.interp(wavelength[weak_si_trough], [si_max_wave_1, si_max_wave_2], [si_max_1, si_max_2])
	interp_flux[strong_si_trough] = np.interp(wavelength[strong_si_trough], [si_max_wave_2, si_max_wave_3], [si_max_2, si_max_3])

	v_strong, si_min_wave = measure_velocity(wavelength,flux, 5900., 6300.)
	v_weak, si_weak_min_wave = measure_weak_si_velocity(wavelength,flux)

	si_min_index = np.where(wavelength == si_min_wave)
	si_weak_min_index = np.where(wavelength == si_weak_min_wave)

	#Line ratio with fractional depths
	# strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])/interp_flux[si_min_index]
	# weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])/interp_flux[si_weak_min_index]

	strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])
	weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])

	ratio = weak_line/strong_line
	ratio = ratio[0]

	# plt.plot(wavelength,flux)
	# plt.plot(wavelength,sm_flux)
	# plt.plot(wavelength, interp_flux)
	# plt.plot(wavelength[si_min_index], interp_flux[si_min_index], 'o', color='orange')
	# plt.plot(wavelength[si_min_index], sm_flux[si_min_index], 'o', color='orange')
	# plt.plot(wavelength[si_weak_min_index], interp_flux[si_weak_min_index], 'o', color='orange')
	# plt.plot(wavelength[si_weak_min_index], sm_flux[si_weak_min_index], 'o', color='orange')
	# plt.show()

	return ratio

def measure_ca_ratio(wavelength,flux, varflux = None, wave1 = 3550., wave2=3680.):
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)

	sm_flux = df.gsmooth(wavelength, flux, varflux, .001)

	# ca_range_1 = np.where((wavelength > 3550.) & (wavelength < 3680.))
	ca_range_1 = np.where((wavelength > wave1) & (wavelength < wave2))
	ca_range_2 = np.where((wavelength > 3890.) & (wavelength < 3970.))

	ca_wave_1 = wavelength[ca_range_1]
	ca_wave_2 = wavelength[ca_range_2]

	ca_flux_1 = sm_flux[ca_range_1]
	ca_flux_2 = sm_flux[ca_range_2]

	ca_max_1 = np.amax(ca_flux_1)
	ca_max_2 = np.amax(ca_flux_2)

	ca_max_index_1 = np.where(sm_flux == ca_max_1)
	ca_max_index_2 = np.where(sm_flux == ca_max_2)

	ca_max_wave_1 = wavelength[ca_max_index_1][0]
	ca_max_wave_2 = wavelength[ca_max_index_2][0]

	print ca_max_2/ca_max_1

	plt.plot(wavelength,flux)
	plt.plot(wavelength,sm_flux)
	plt.plot(wavelength[ca_max_index_1], sm_flux[ca_max_index_1], 'o', color='orange')
	plt.plot(wavelength[ca_max_index_2], sm_flux[ca_max_index_2], 'o', color='orange')
	plt.show()

	return ca_max_2/ca_max_1

def measure_velocity(wavelength, flux, wave1, wave2, vexp=.001, rest_wave=6355., varflux=None):
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)

	sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
	si_range = np.where((wavelength > wave1) & (wavelength < wave2))
	si_wave = wavelength[si_range]
	if len(si_wave) == 0:
		return np.nan, np.nan
	si_flux = sm_flux[si_range]
	si_min = np.amin(si_flux)
	si_min_index = np.where(si_flux == si_min)
	si_min_wave = si_wave[si_min_index][0]

	c = 299792. # km/s
	# rest_wave = 6355. #Angstroms

	v = c*((rest_wave/si_min_wave)**2. - 1)/(1+((rest_wave/si_min_wave)**2.))
	# plt.plot(wavelength, flux)
	# plt.plot(wavelength, sm_flux)
	# plt.plot(si_min_wave, si_min, 'o', color='orange')
	# plt.show()
	return v, si_min_wave

def measure_weak_si_velocity(wavelength, flux, varflux=None):
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)

	varflux = np.zeros(len(wavelength), float)
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
	return v, si_min_wave

	#carbon II is 6580A (2006bt good example)

def measure_si_velocity_from_raw(wavelength,flux,z, varflux=None):
	# varflux = 1./SN.ivar[SN.x1:SN.x2]
	# plt.plot(wavelength,flux)
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)

	wavelength = wavelength/(1.+z)
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

def measure_C_velocity_from_raw(wavelength,flux,z, varflux=None):
	# varflux = 1./SN.ivar[SN.x1:SN.x2]
	# plt.plot(wavelength,flux)
	if varflux == None:
		varflux = np.zeros(len(wavelength), float)
	wavelength = wavelength/(1.+z)
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

def resample_spectrum(SN):
	wavelength = SN.wavelength[SN.x1:SN.x2]
	flux = SN.flux[SN.x1:SN.x2]
	varflux = 1./SN.ivar[SN.x1:SN.x2]

# spec_file = r'..\..\SNIa_maximum_light.flm'
# with open(spec_file) as spec:
# 	spectrum = np.transpose(np.loadtxt(spec))
# 	wavelength = spectrum[0]
# 	flux = spectrum[1]
# 	measure_si_ratio(wavelength,flux)
# 	measure_ca_ratio(wavelength,flux)
	# z = 0.005274
	# # print  measure_si_velocity_from_raw(wavelength, flux, z)
	# print  measure_C_velocity_from_raw(wavelength, flux, z)