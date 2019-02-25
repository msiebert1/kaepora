import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import magnitudes as mag
import find_event_data as fed
import magnitudes as mag
import matplotlib.pyplot as plt
import photometry as ph

import prep
from astropy import units as u
from specutils import Spectrum1D
import test_dered
import copy
import scipy.signal as sig

# np.set_printoptions(threshold=np.nan)

SN_Array = []
full_array = []
compare_spectrum = []

def reject_outliers(data, m=1):
	for i in range(len(data)):
		if abs(data[i] - np.mean(data)) > m * np.std(data):
			data[i] = 0.
	
	return data

def specific_query(test_query):
	# event = "'" + event + "'"
	query = "SELECT * FROM Supernovae " + test_query
	print query
	SN_Array = fed.grab_event_data(query)
	return SN_Array

def find_missing_data(SN):

	# sm_flux = sig.savgol_filter(SN.flux, 3, 1)
	# plt.plot(SN.wavelength, sm_flux)
	# plt.plot(SN.wavelength, SN.flux)
	# plt.show()

	dw = SN.wavelength[1] - SN.wavelength[0]
	df = np.diff(SN.flux)
	df_dw = df/dw

	dw2 = dw**2
	df2 = np.diff(SN.flux, 2)
	df2_dw2 = df2/dw2

	# print df2_dw2[np.where((SN.wavelength >= 9000.) & (SN.wavelength <= 9090.))]
	tol = 1.e-5
	avg_flux_mag = np.average(np.absolute(SN.flux[SN.x1:SN.x2]))

	mis_data = []
	old_ivar = copy.copy(SN.ivar)
	# for i in range(len(df2_dw2)):
	# 	if np.absolute(df2_dw2[i])/avg_flux_mag < tol and i > SN.x1 and i < SN.x2:
	# 		mis_data.append(i)
	# 		SN.ivar[i] = None

	non_zero_ivar = np.where(SN.ivar != 0.)
	avg_ivar = np.average(SN.ivar[non_zero_ivar])
	sig_ivar = np.std(SN.ivar[non_zero_ivar])

	for i in range(len(SN.ivar)):
		if SN.ivar[i] > 2.*sig_ivar:
			mis_data.append(i)
			SN.ivar[i-20:i+20] = None


	# SN.ivar[mis_data] = None
	print mis_data
	print np.average(np.absolute(SN.flux[SN.x1:SN.x2]))

	# print min(np.absolute(df2_dw2[SN.x1:SN.x2]))
	# print SN.wavelength[np.where(np.absolute(df2_dw2) == min(np.absolute(df2_dw2[SN.x1:SN.x2])))]

	plt.subplot(2,1,1)
	# plt.plot(SN.wavelength[SN.x1:SN.x2], sm_flux[SN.x1:SN.x2])
	plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
	plt.subplot(2,1,2)
	plt.plot(SN.wavelength[SN.x1:SN.x2], old_ivar[SN.x1:SN.x2], 'r')
	# plt.plot(SN.wavelength[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2])
	# plt.subplot(4,1,3)
	# plt.plot(SN.wavelength[SN.x1:SN.x2], df_dw[SN.x1:SN.x2])
	# plt.subplot(4,1,4)
	# plt.plot(SN.wavelength[SN.x1:SN.x2], df2_dw2[SN.x1:SN.x2])
	plt.show()



# SN_Array = specific_query("where filename = 'sn1994t-19940715-ui.flm'")
# SN_Array = specific_query("where filename = 'sn1994S-19940612.26-mmt.flm'")
# SN_Array = specific_query("where phase  between -2 and 2 and dm15 < 1.1 and source != 'cfa'")

# for SN in SN_Array:
# 	zero_data = np.where(SN.ivar == 0)
# 	zero_data = np.array(zero_data[0])
# 	sub_lists = get_sub_list(zero_data)
# 	sub_lists = sub_lists[1:-1]
# 	for l in sub_lists:
# 		fill_val = (SN.ivar[l[0]-1] + SN.ivar[l[-1]+1])/2.
# 		for i in l:
# 			SN.ivar[i] = fill_val

# SN = SN_Array[0]

# find_missing_data(SN)

# for SN in SN_Array:
# 	pre_scale = (1.e-15/np.average(SN.flux[SN.x1:SN.x2]))
# 	SN.flux = pre_scale*SN.flux
# 	SN.ivar = SN.ivar/(pre_scale*pre_scale)
# 	host_reddened = prep.ReadExtin('../data/info_files/ryan_av.txt')
# 	old_wave = SN.wavelength*u.Angstrom        # wavelengths
# 	old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
# 	old_ivar = SN.ivar*u.Unit('W m-2 angstrom-1 sr-1')
# 	spec1d = Spectrum1D.from_array(old_wave, old_flux)
# 	new_flux, new_ivar, corrected = test_dered.host_correction(host_reddened, SN.name, old_wave, old_flux, old_ivar)
# 	SN.flux = new_flux.value
# 	SN.ivar = new_ivar.value

# 	if corrected:
# 		print SN.filename, SN.source
# 		plt.subplot(2,1,1)
# 		plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
# 		plt.subplot(2,1,2)
# 		plt.plot(SN.wavelength[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2])
# 		# plt.plot(SN.wavelength[SN.x1:SN.x2], 1./(SN.ivar[SN.x1:SN.x2])**.5)
# 		plt.show()


# for SN in SN_Array:
# 	print SN.name, SN.source
# SN_Array = specific_query("where SN = '2003du'")
# for b in SN_Array[0].phot:
# 	print b
# for SN in SN_Array:
# 	print SN.filename, SN.phase

# phot = []
# scales = []
# for SN in SN_Array:
# 	valid_bands = mag.valid_bands(SN)
# 	print valid_bands
# 	scale, mags_from_phot = mag.scale_flux_to_photometry(SN, valid_bands)
# 	scales.append(scale)
# 	phot.append(mags_from_phot)
# 	SN.flux *= scale
# 	plt.plot(SN.wavelength, SN.flux)
# 	plt.show()

##Plots all photometry for a given event
p = SN_Array[0].phot
times = []
ms = []
for band in p:
	t = []
	m = []
	for i in range(len(p[band][0])):
		t.append(float(p[band][0][i]))
		m.append(float(p[band][1][i][0]))
	times.append(t)
	ms.append(m)
for i in range(len(times)):
	plt.plot(times[i], ms[i], 'o')
plt.gca().invert_yaxis()
plt.show()

for i in range(len(SN_Array)):
	print SN_Array[i].phot
	phot = SN_Array[i].phot
	if phot[i] is not None and len(phot[i]) >= 2:
		if 'R' in phot[i] and 'V' in phot[i] and phot[i]['R'] is not None and phot[i]['V'] is not None:
			color_from_phot = phot[i]['V'] - phot[i]['R']
			color_from_spec = mag.ab_mag(SN_Array[i], 'V') - mag.ab_mag(SN_Array[i], 'R')
			print scales[i], color_from_spec, color_from_phot


valid_bands = mag.valid_bands(SN_Array[10])
mag.scale_flux_to_photometry(SN_Array[10], valid_bands)