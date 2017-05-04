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

def specific_query(test_query):
	# event = "'" + event + "'"
	query = "SELECT * FROM Supernovae " + test_query
	print query
	SN_Array = fed.grab_event_data(query)
	return SN_Array

# SN_Array = specific_query("where filename = 'sn1994t-19940715-ui.flm'")
SN_Array = specific_query("where filename = 'sn1994S-19940612.26-mmt.flm'")
# SN_Array = specific_query("where phase  between -2 and 2 and dm15 < 1.1 and source != 'cfa'")

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


# for i in range(len(SN_Array)):
# 	phot = SN_Array[i].phot
# 	if 'R' in phot and 'V' in phot and phot['R'] is not None and phot['V'] is not None:
# 		color_from_phot = phot['V'] - phot['R']
# 		color_from_spec = mag.ab_mag(SN_Array[i], 'V') - mag.ab_mag(SN_Array[i], 'R')
# 		print scales[i], color_from_spec, color_from_phot


valid_bands = mag.valid_bands(SN_Array[0])
mag.scale_flux_to_photometry(SN_Array[0], valid_bands)