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

questionable_files = []
# 1. Read in data from data_notes.txt, check out python's readlines() function
# 2. Use loops and if statements to store all questionable filenames in the above list
# (the list should replace the "files" list below. Note the extra quotation marks.)

files = ["'SN08bf_080328_b01_NTT_EM.dat'", "'SN05ke_051123_b01_DUP_MS.dat'"] # example file array 
f_string = "("
for f in files:
	f_string = f_string + f + ','
f_string = f_string[0:-1] + ")"

print f_string
SN_Array = specific_query("where filename in "+ f_string)
for SN in SN_Array:
	plt.plot(SN.wavelength, SN.flux)
	plt.show()