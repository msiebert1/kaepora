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
							all_phot[band][0].append(float(phot.get(time_key, None)))
							all_phot[band][1].append(float(phot.get(mag_key, None)))

		if len(all_phot) == 0:
			all_phot = None
		return all_phot	
	else:
		return None	

def find_event_in_osc(sn_name):
	files = glob.glob("..\..\osc_data\osc_overlap\*.json")
	for f in files:
		if sn_name + '.json' in f.lower():
			return f
	return None

