import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import magnitudes as mag
import find_event_data as fed
import magnitudes as mag
import matplotlib.pyplot as plt

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

SN_Array = specific_query("where phase  between -1 and 1")

for SN in SN_Array:
	plt.plot(SN.wavelength, SN.flux)
	plt.show()



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
