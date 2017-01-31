import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import magnitudes as mag
import find_event_data as fed
import magnitudes as mag
import matplotlib.pyplot as plt
import photometry as ph

# np.set_printoptions(threshold=np.nan)

SN_Array = []
full_array = []
compare_spectrum = []

def specific_query(test_query):
	# event = "'" + event + "'"
	query = "SELECT * FROM Supernovae " + test_query
	print query
	SN_Array = fed.grab_event_data(query)
	return SN_Array

# SN_Array = specific_query("where SN = '2002cx'")
# phot = []
# scales = []
# print SN_Array[0].phot
all_phot = ph.get_photometry('2002cx')
print all_phot
# for SN in SN_Array:
# 	valid_bands = mag.valid_bands(SN)
# 	scale, mags_from_phot = mag.scale_flux_to_photometry(SN, valid_bands)
# 	scales.append(scale)
# 	phot.append(mags_from_phot)
# 	SN.flux *= scale

# for i in range(len(SN_Array)):
# 	if phot[i] is not None and len(phot[i]) >= 2:
# 		if 'R' in phot[i] and 'V' in phot[i] and phot[i]['R'] is not None and phot[i]['V'] is not None:
# 			color_from_phot = phot[i]['V'] - phot[i]['R']
# 			color_from_spec = mag.ab_mag(SN_Array[i], 'V') - mag.ab_mag(SN_Array[i], 'R')
# 			print scales[i], color_from_spec, color_from_phot


# valid_bands = mag.valid_bands(SN_Array[10])
# mag.scale_flux_to_photometry(SN_Array[10], valid_bands)



# SN_Array = specific_query("source = 'bsnip'")
# for SN in SN_Array:
# 	print SN.mjd


#Photometry stuff
# SN_Array = specific_query("SN = '2005cf'")
# SN = SN_Array[30]
# valid_bands = ["U", "B", "V", "R", "I"]
# mag.scale_flux_to_valid_bands(SN, valid_bands)




# scale = mag.scale_flux_to_Vband(SN)


# times = np.asarray(vband[0])
# times = times.astype(np.float)
# mags = []
# for datum in vband[1]:
# 	mags.append(float(datum[0]))

# plt.plot(times, mags)
# plt.gca().invert_yaxis()
# plt.show()
