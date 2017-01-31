import numpy as np
import json
from pprint import pprint
import glob
import matplotlib.pyplot as plt

def get_photometry(sn_name):

	file = find_event_in_osc(sn_name)

	if file != None:
		phot_key = "photometry"
		band_key = "band"
		time_key = "time"

		mag_key = "magnitude"
		e_key = "e_magnitude"
		e_low_key = "e_lower_magnitude"
		e_up_key = "e_upper_magnitude"
		zero_point_key = "zeropoint"
		sys_key = "system"
		upper_limit_key = "upperlimit"
		limit_sigma_key = "upperlimitsigma"
		kcorr_key = "kcorrected"
		hcorr_key = "scorrected"
		mcorr_key = "mcorrected"

		all_phot = {}

		with open (file) as f:
			data = json.load(f)
			filename = file
			file_key = file.split('\\')[-1][:-5]
			data = data.get(file_key)
			if phot_key in data:
				for phot in data[phot_key]:
					if phot.get(band_key, None) != None:
						if phot.get(band_key) not in all_phot:
							all_phot[phot.get(band_key, None)] = [[],[]]
						band = phot.get(band_key, None)
						if phot.get(time_key, None) != None and phot.get(mag_key, None) != None:
							all_phot[band][0].append(phot.get(time_key, None))
							#Each of the following fields is populated if the data exist (None if no data).
							#currently: 0 - magnitude, 1 - error, 2 - lower error, 3 - upper error, 4 - zero point, 5 - photometric system,
							#			6 - is upper limit, 7 - upper limit sigma, 8 - k corrected, 9 - host corrected, 10 - MW corrected
							# all_phot[band][1].append([phot.get(mag_key, None), phot.get(e_key, None), phot.get(e_low_key, None),
							# 						  phot.get(e_up_key, None), phot.get(zero_point_key, None), phot.get(sys_key, None),
							# 						  phot.get(upper_limit_key, None), phot.get(limit_sigma_key, None), phot.get(kcorr_key, None),
							# 						  phot.get(hcorr_key,None), phot.get(mcorr_key, None)])
							all_phot[band][1].append([phot.get(mag_key, None), phot])
		if len(all_phot) == 0:
			all_phot = None
		return all_phot	
	else:
		return None	

def find_event_in_osc(sn_name):
	files = glob.glob("..\..\osc_data\osc_overlap\*.json")
	for f in files:
		if sn_name.lower() + '.json' in f.lower():
			return f
	return None

