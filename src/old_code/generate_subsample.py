import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import magnitudes as mag
import find_event_data as fed
import magnitudes as mag
import matplotlib.pyplot as plt
import jdcal

# np.set_printoptions(threshold=np.nan)

SN_Array = []
full_array = []
compare_spectrum = []

def specific_query(test_query):
	# event = "'" + event + "'"
	query = "SELECT * FROM Supernovae where " + test_query
	print query
	SN_Array = fed.grab_event_data(query)
	return SN_Array

# SN_Array = specific_query("Phase between -4 and 4 and dm15 >= 0.0")
SN_Array = specific_query("Phase between -4 and 4 and dm15 >= 0.0")
new_arr = []
short_specs = []
for SN in SN_Array:
	if SN.wavelength[SN.x1] <= 4000. and SN.wavelength[SN.x2] >= 7000.:
		new_arr.append(SN)
		if SN.wavelength[SN.x2] >= 9000.:
			short_specs.append(SN)


print len(new_arr), len(short_specs)


files_to_remove = ['sn1996ai-19960623.27-fast.flm','sn1997dt-19971204.11-fast.flm',
				   'sn1997dt-19971206.11-fast.flm', 'sn2003cg-20030329.27-mmt.flm',
				   'sn2003ch-20030330.16-fast.flm', 'sn2003gq-20030729.39-fast.flm', 
				   'sn2004ef-20040915.30-fast.flm', 'sn2005A-20050107.25-fast.flm', 
				   'sn2005A-20050108.13-fast.flm', 'sn2005bo-20050418.21-fast.flm',
				   'sn2005cc-20050531.20-fast.flm', 'sn2005cc-20050601.23-fast.flm',
				   'sn2005eq-20051009.42-fast.flm', 'sn2006az-20060403.30-fast.flm',
				   'sn2006br-20060427.33-fast.flm', 'sn2006cm-20060528.45-fast.flm',
				   'sn2007S-20070209.30-fast.flm', 'sn2007S-20070210.32-fast.flm',
				   'sn2007S-20070212.33-fast.flm', 'sn2008ar-20080307.38-fast.flm',
				   'sn2000cp-20000627-ui-corrected.flm', 'sn2000dn-20001006-uri-corrected.flm',
				   'sn2006cj-20060528.349-ui.flm', 'SN07S_070212_b01_DUP_WF.dat',
				   '2005cf_20050601_3243_9720_00.dat', 'sn1996ai-19960620-uo.flm',
				   'sn1996ai-19960621-uo.flm', 'sn2003gn-20030805.271-ui.flm',
				   'sn2003gq-20030727-ui.flm', 'SN04gs_041213_b01_DUP_WF.dat',
				   'SN08bf_080326_b01_CLA_LD.dat', '2002bo_20020320_3305_9125_00.dat',
				   '2002bo_20020321_3306_9127_00.dat']

# plt.figure(num = 2, dpi = 100, figsize = [20, 10], facecolor = 'w')
# i = 1
# for SN in new_arr:
# 	if SN.filename not in files_to_remove:
# 		print i, SN.filename, SN.phase, SN.dm15
# 		plt.figure(num = 2, dpi = 100, figsize = [20, 10], facecolor = 'w')
# 		plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
# 		plt.show()
# 		i +=1

final_sample = []
for SN in new_arr:
	if SN.filename not in files_to_remove:
		final_sample.append(SN)

final_short_sample = []
for SN in short_specs:
	if SN.filename not in files_to_remove:
		final_short_sample.append(SN)

print len(final_sample), len(final_short_sample)

# for SN in final_sample:
# 	if SN.phase >= 0.:
# 		sign = 'p'
# 	else:
# 		sign = 'm'
# 	dm15_str = str(SN.dm15)
# 	abs_phase = np.absolute(SN.phase)
# 	phase_str = str(abs_phase)
# 	with open('../../Sample_for_Ryan/' + 'sn' + SN.name + '-' + sign + phase_str + '-' + dm15_str + '.flm', 'w') as file:
# 		wave = np.array(SN.wavelength[SN.x1:SN.x2])
# 		flux = np.array(SN.flux[SN.x1:SN.x2])
# 		data = np.array([wave,flux])
# 		data = data.T
# 		np.savetxt(file, data)

for SN in final_short_sample:
	if SN.phase >= 0.:
		sign = 'p'
	else:
		sign = 'm'
	dm15_str = str(SN.dm15)
	abs_phase = np.absolute(SN.phase)
	phase_str = str(abs_phase)
	with open('../../Sample_for_Ryan/' + 'sn' + SN.name + '-' + sign + phase_str + '-' + dm15_str + '.flm', 'w') as file:
		wave = np.array(SN.wavelength[SN.x1:SN.x2])
		flux = np.array(SN.flux[SN.x1:SN.x2])
		data = np.array([wave,flux])
		data = data.T
		np.savetxt(file, data)


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
