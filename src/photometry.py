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

def get_csp_photometry(sn_name):
	file = find_csp_phot(sn_name)
	phot_dict = {}
	index_dict = {}
	if file != None:
		with open(file) as f:
			lines = f.readlines()
			for i in range(len(lines[4].split()[1:])):
				band = lines[4].split()[1:][i]
				if band != "+/-" and band != "MJD":
					phot_dict[band] = [[],[],[]]
					index_dict[band] = i
			for line in lines[5:]:
				for b in phot_dict.keys():
					if line.split()[index_dict[b]] != '99.900':
						phot_dict[b][0].append(float(line.split()[0]))
						phot_dict[b][1].append(float(line.split()[index_dict[b]]))
						phot_dict[b][2].append(float(line.split()[index_dict[b] + 1]))
			return phot_dict
	else:
		return None

def find_csp_phot(sn_name):
	files = glob.glob("..\..\csp_photometry\CSP_Photometry_DR2\*.dat")
	for f in files:
		if sn_name.lower() in f.lower():
			return f
	return None

def find_event_in_osc(sn_name):
	files = glob.glob("..\..\osc_data\osc_sns_new\*.json")
	for f in files:
		if sn_name.lower() + '.json' in f.lower():
			return f
	return None