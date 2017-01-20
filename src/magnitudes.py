import matplotlib.pyplot as plt
import numpy as np
import speclite.filters as filters
import astropy.units as u
import scipy.optimize as opt

def ab_mags(SN):
	phases = []
	mags = []
	# bessell = filters.load_filters('bessell-*')
	Uband = filters.load_filter('bessell-U')
	Bband = filters.load_filter('bessell-B')
	Vband = filters.load_filter('bessell-V')
	Rband = filters.load_filter('bessell-R')
	Iband = filters.load_filter('bessell-I')

	wave = SN.wavelength
	flux = SN.flux

	# phases.append(SN.phase)
	Umag = Uband.get_ab_magnitude(flux, wave)
	Bmag = Bband.get_ab_magnitude(flux, wave)
	Vmag = Vband.get_ab_magnitude(flux, wave)
	Rmag = Rband.get_ab_magnitude(flux, wave)
	Imag = Iband.get_ab_magnitude(flux, wave)
	mags.append([[SN.phase],[Umag, Bmag, Vmag, Rmag, Imag]])

	return mags

def valid_bands(SN):
	U_range = [3000.,4500.]
	B_range = [3500.,5800.]
	V_range = [4500.,7000.]
	R_range = [5400.,9000.]
	I_range = [7000.,9500.]
	min_wave = SN.wavelength[SN.x1]
	max_wave = SN.wavelength[SN.x2]


def scale_flux_to_valid_bands(SN, valid_bands):
	filts = []
	for b in valid_bands:
		filts.append(filters.load_filter('bessell-' + b))

	if filts[0].get_ab_magnitude(SN.flux, SN.wavelength) < 0.:
		SN.flux = SN.flux*1.e-15

	guess = 1.
	mags_from_phot = []
	mag = 0.
	for i in range(len(valid_bands)):
		phot = SN.phot.get(valid_bands[i], None)
		if phot is not None:
			min_diff = np.absolute(SN.mjd - float(phot[0][0]))
			###FIX date diffs too big!!!
			print min_diff
			for j in range(len(phot[0])):
				diff = np.absolute(SN.mjd - float(phot[0][i]))
				if diff < min_diff:
					min_diff = diff
					mag = float(phot[1][i][0])
			mags_from_phot.append(mag)
		else:
			mags_from_phot.append(None)

	scale = opt.minimize(total_phot_offset, guess, args = (SN, filts, mags_from_phot)).x

	final = filts[0].get_ab_magnitude(scale*SN.flux, SN.wavelength)
	print scale[0], final, mags_from_phot[0]


def total_phot_offset(scale, SN, filters, target_mags):
	diff_sum = 0.
	for i in range(len(filters)):
		if target_mags[i] is not None:
			mag_from_spec = filters[i].get_ab_magnitude(scale*SN.flux, SN.wavelength)
			diff = (mag_from_spec - target_mags[i])**2.
			diff_sum += diff
	return diff_sum


def scale_flux_to_Vband(SN):
	Vband = filters.load_filter('bessell-V')

	Vband_phot = SN.phot.get("V", None)

	if Vband_phot is None:
		return None
	if Vband.get_ab_magnitude(SN.flux, SN.wavelength) < 0.:
		print "flux has been scaled"
		SN.flux = SN.flux*1.e-15

	guess = 1.
	min_diff = np.absolute(SN.mjd - float(Vband_phot[0][0]))
	v_mag = 0.
	for i in range(len(Vband_phot[0])):
		diff = np.absolute(SN.mjd - float(Vband_phot[0][i]))
		if diff < min_diff:
			min_diff = diff
			Vmag_from_phot = float(Vband_phot[1][i][0])

	scale = opt.minimize(phot_offset_using_flux, guess, args = (SN, Vband, Vmag_from_phot)).x

	final = Vband.get_ab_magnitude(scale*SN.flux, SN.wavelength)
	print scale[0], final, Vmag_from_phot

	return scale

def phot_offset_using_flux(scale, SN, band, target_mag):
	if scale < 0.:
		return 1.e6
	mag_from_spec = band.get_ab_magnitude(scale*SN.flux, SN.wavelength)
	return np.absolute(mag_from_spec - target_mag)

def generate_event_list(SN_Array):
	unique_events = []
	event_file = open("event_list.txt", 'w')
	for SN in SN_Array:
		if SN.name not in unique_events:
			unique_events.append(SN.name)
			event_file.write(SN.name + ' | '+ '\n')
	event_file.close()

