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
f = open('data_notes_csp_other.txt', 'r')
questionable_files = []
for line in f.readlines():
    comment = line.split()[2]
    if 't' in comment.split(','):
        questionable_files.append("'" + line.split()[1] + "'")


# questionable_files = ["'sn1989m-19900501-opt.flm'", 	#interpolated
# 				  "'sn1990y-19900830-opt2.flm'",	#strong emission
# 				  "'sn1991bg-19920523-uoi.flm'",	#contaminated
# 				  "'sn1991s-19910505-final.flm'",	#very low SNR
# 				  "'sn1991t-19920411-uoi.flm'", 	#contaminated
# 				  "'sn1992m-19920313-uv.flm'", 		#low SNR, contaminated
# 				  "'sn1994t-19940715-ui.flm'", 		#oscillation, systematic?
# 				  "'sn2000dn-20001006-uri-corrected.flm'", #telluric
# 				  "'sn2002bf-20020311-ui-corrected.flm'", #noisy, lines systematic?
# 				  "'sn2002cv-20020608-ui.flm'", #red and blue side different SNRs
# 				  "'sn2004bv-20040711-ui-corrected.flm'", #telluric
# 				  "'sn2005gj-20051202.403-deimos.flm'" #strong emission
# 				  "'sn1994S-19940612.26-mmt.flm'"]	#silicon line interpolated


# questionable_files = ["'sn1992g-19920625-final.flm'",
#                       "'sn1992m-19920313-ir.flm'",
#                       "'sn1992m-19920313-uv.flm'",
#                       "'sn1993ab-19931022-ui.flm'",
#                       "'sn1993ac-19931022-uoi.flm'",
#                       "'sn1993ac-19931108-ui.flm'",
#                       "'sn1993ae-19931117-uoi.flm'",
#                       "'sn1994d-19940714.flm'",
#                       "'sn1994e-19940313-opt.flm'",
#                       "'sn1994j-19940404-opt.flm'",
#                       "'sn1994q-19940714-uoi.flm'",
#                       "'sn1994q-19940804-ui.flm'",
#                       "'sn1994t-19940715-ui.flm'",
#                       "'sn1994x-19940901-ui.flm'",
#                       "'sn1995a-19950224-opt-blo.flm'",
#                       "'sn1995ac-19950926-uoi.flm'",
#                       "'sn2000cx-20011022.270-joined.flm'",
#                       "'sn2000dd-20000827-ui.flm'",
#                       "'sn2000df-20000826-ui.flm'",
#                       "'sn2004bv-20040711-ui-corrected.flm'",
#                       "'sn2006en-20060919.209-ui-corrected.flm'",
#                       "'sn2006en-20060925.232-ui-corrected.flm'",
#                       "'sn2006es-20060902.482-ui.flm'",
#                       "'sn2006es-20060919.499-ui.flm'",
#                       "'sn2006ev-20060919.183-ui.flm'",
#                       "'sn2006ev-20060925.175-ui.flm'",
#                       "'sn2006gj-20060925.482-ui.flm'",
#                       "'sn2006gr-20060925.202-ui.flm'",
#                       "'sn2006gr-20061024.150-ui.flm'",
#                       "'sn2006gt-20060925.348-ui.flm'",
#                       "'sn2006ha-20061124.351-600b.flm'",
#                       "'sn2006hb-20061030.489-ui.flm'",
#                       "'sn2006je-20061024.247-ui-correct'",
#                       "'sn2006je-20061030.313-ui.flm'",
#                       "'sn2006ke-20061024.397-ui.flm'",
#                       "'sn1996X-19960419.30-fast.flm'",
#                       "'sn1996X-19960420.31-fast.flm'",
#                       "'sn1996X-19960426.28-fast.flm'",
#                       "'sn1996Z-19960522.17-fast.flm'",
#                       "'sn1996Z-19960524.18-fast.flm'",
#                       "'sn1996Z-19960528.16-fast.flm'"]

# questionable_files = ["'1996X_19960614_2806_9811_00.dat'"]
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