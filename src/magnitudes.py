import matplotlib.pyplot as plt
import numpy as np
import speclite.filters as filters
import astropy.units as u
import scipy.optimize as opt
from scipy.interpolate import interp1d

def ab_mag(SN, filt):
	f = filters.load_filter('bessell-' + filt)
	mag = f.get_ab_magnitude(SN.flux, SN.wavelength)
	return mag


def valid_bands(SN):
	valid_bands = []
	band_waves = {"U": [3000.,4500.], 
				  "B": [3500.,5800.], 
				  "V": [4500.,7000.], 
				  "R": [5400.,9000.], 
				  "I": [7000.,9500.]}
	min_wave = SN.wavelength[SN.x1]
	max_wave = SN.wavelength[SN.x2]

	for b in band_waves:
		if min_wave < band_waves[b][0] and max_wave > band_waves[b][1]:
			valid_bands.append(b)

	return valid_bands


def interp_LC(phot):
	min_points = 5 #what should this be?
	temp = np.asarray(phot[0])
	times = temp.astype(float)
	if len(times) > min_points:
		mags = []
		for pt in phot[1]:
			mags.append(float(pt[0]))

		lc = interp1d(times, mags, bounds_error = True)
		return lc
	else:
		return None


def generate_photometry_for_epoch(SN, valid_bands):
	all_phot = SN.phot
	light_curves = {}
	for band in valid_bands:
		if band in all_phot:
			light_curves[band] = interp_LC(all_phot[band])
		else:
			light_curves[band] = None

	mags_from_phot = {}
	for band in light_curves:
		if light_curves[band] is not None:
			try:
				mags_from_phot[band] = light_curves[band](SN.mjd) #check for extrap. this will throw error
			except ValueError:
				mags_from_phot[band] = None
		else:
			mags_from_phot[band] = None

	return mags_from_phot


def scale_flux_to_photometry(SN, valid_bands):

	if len(valid_bands) > 0:
		mags_from_phot = generate_photometry_for_epoch(SN, valid_bands)

		filts = {}
		for b in valid_bands:
			filts[b] = filters.load_filter('bessell-' + b)

		if filts[valid_bands[0]].get_ab_magnitude(SN.flux, SN.wavelength) < 0.: #temporary fix for previously scaled fluxes
			SN.flux = SN.flux*1.e-15

		guess = 1.
		scale = opt.minimize(total_phot_offset, guess, args = (SN, filts, mags_from_phot)).x
		scale = scale[0]
		final = filts[valid_bands[0]].get_ab_magnitude(scale*SN.flux, SN.wavelength)
		# print scale, final, mags_from_phot #for testing output
	else:
		scale = 1.
		mags_from_phot = None

	return scale, mags_from_phot


def total_phot_offset(scale, SN, filters, target_mags):
	diff_sum = 0.
	if scale < 0.: #do not allow negative scale values
		return 1.e6
	else:
		for f in filters:
			if target_mags[f] is not None:
				mag_from_spec = filters[f].get_ab_magnitude(scale*SN.flux, SN.wavelength)
				diff = (mag_from_spec - target_mags[f])**2.
				diff_sum += diff
		return diff_sum


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

