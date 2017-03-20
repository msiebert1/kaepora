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
f = open('data_notes.txt', 'r')
questionable_files = []
for line in f.readlines():
	comment = line.split()[2]

	if comment.split(',')[0] is 'q':
		questionable_files.append("'" + line.split()[1] + "'")

questionable_files = ["'sn1989m-19900501-opt.flm'", 	#interpolated
				  "'sn1990y-19900830-opt2.flm'",	#strong emission
				  "'sn1991bg-19920523-uoi.flm'",	#contaminated
				  "'sn1991s-19910505-final.flm'",	#very low SNR
				  "'sn1991t-19920411-uoi.flm'", 	#contaminated
				  "'sn1992m-19920313-uv.flm'", 		#low SNR, contaminated
				  "'sn1994t-19940715-ui.flm'", 		#oscillation, systematic?
				  "'sn2000dn-20001006-uri-corrected.flm'", #telluric
				  "'sn2002bf-20020311-ui-corrected.flm'", #noisy, lines systematic?
				  "'sn2002cv-20020608-ui.flm'", #red and blue side different SNRs
				  "'sn2004bv-20040711-ui-corrected.flm'", #telluric
				  "'sn2005gj-20051202.403-deimos.flm'" #strong emission
				  "'sn1994S-19940612.26-mmt.flm'"]	#silicon line interpolated

f_string = "("
for f in questionable_files:
	f_string = f_string + f + ','
f_string = f_string[0:-1] + ")"

print f_string

SN_Array = specific_query("where filename in "+ f_string)
for SN in SN_Array:
	print SN.filename, SN.phase, SN.SNR
	plt.figure(num = 2, dpi = 100, figsize = [30, 15], facecolor = 'w')
	plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
	plt.show()

bad_files = ['sn1991t-19910820-ir2.flm',	#interpolated
			 'sn1989m-19900501-opt.flm',	#interpolated
			 'sn1991t-19920313-ui.flm',		#interpolated
			 'sn1992m-19920313-ir.flm', 	#interpolated
			 'sn1994e-19940313-opt.flm', 	#interpolated
			 'sn2000cx-20011022.270-joined.flm', #interpolated
			 ]