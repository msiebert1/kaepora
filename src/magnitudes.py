import matplotlib.pyplot as plt
import numpy as np
import speclite.filters as filters
import astropy.units as u

def ab_mags(SN_Array):
	phases = []
	mags = []
	sdss = filters.load_filters('sdss2010-*')
	for SN in SN_Array:
		wave = SN.wavelength
		flux = SN.flux
		# phases.append(SN.phase)
		mag_table = sdss.get_ab_magnitudes(flux, wave)
		mags.append([[SN.phase],[mag_table['sdss2010-u'].quantity[0], mag_table['sdss2010-g'].quantity[0], 
					 mag_table['sdss2010-r'].quantity[0], mag_table['sdss2010-i'].quantity[0],
					 mag_table['sdss2010-z'].quantity[0]]])

	return mags

def generate_event_list(SN_Array):
	unique_events = []
	event_file = open("event_list.txt", 'w')
	for SN in SN_Array:
		if SN.name not in unique_events:
			unique_events.append(SN.name)
			event_file.write(SN.name + ' | '+ '\n')
	event_file.close()
