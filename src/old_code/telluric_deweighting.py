import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import find_event_data as fed
import matplotlib.pyplot as plt

def specific_query(test_query):
	# event = "'" + event + "'"
	query = "SELECT * FROM Supernovae " + test_query
	print query
	SN_Array = fed.grab_event_data(query)
	return SN_Array

def deweight_tell_ranges(SN_Array):
	A_band = [5000,6000]
	B_band = [8000,9000]
	interp_size = 100
	weight_val = 0.
	for SN in SN_Array:
		A_range = np.where((SN.wavelength > A_band[0]) & (SN.wavelength < A_band[1]))
		B_range = np.where((SN.wavelength > B_band[0]) & (SN.wavelength < B_band[1]))
		SN.ivar[A_range] = weight_val
		SN.ivar[B_range] = weight_val

		#continue here

		
	return SN_Array


questionable_files = []

# 1. Read in data from data_notes.txt, check out python's readlines() function
# 2. Use loops and if statements to store all questionable filenames in the above list
# (the list should replace the "files" list below. Note the extra quotation marks.)
f = open('data_notes_csp_other.txt', 'r')
questionable_files = []
for line in f.readlines():
    comment = line.split()[2]
    if 't' in comment.split(','):
        questionable_files.append("'" + line.split()[1] + "'")

f_string = "("
for f in questionable_files:
	f_string = f_string + f + ','
f_string = f_string[0:-1] + ")"

SN_Array = specific_query("where filename in "+ f_string)

deweighted_SN_Array = deweight_tell_ranges(SN_Array)

plt.plot(SN_Array[0].wavelength[SN_Array[0].x1:SN_Array[0].x2], SN_Array[0].ivar[SN_Array[0].x1:SN_Array[0].x2])
plt.show()
