import numpy as np
import os
from astropy.table import Table
from astropy.io import ascii
import scipy.interpolate as inter
import glob
import matplotlib.pyplot as plt
import math

spectra_folders = glob.glob("../../../data/cfa/*")
sn_folders=[]
bad_folders=[]

for i in range(len(spectra_folders)):
	if '.dat' in spectra_folders[i]:
		bad_folders.append(spectra_folders[i])
	elif 'snf' in spectra_folders[i]:
		bad_folders.append(spectra_folders[i])
	else:
		sn_folders.append(spectra_folders[i])

def find_max(file_names):
	parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
	sn_name = parameters["col1"]
	jdate=parameters["col3"]
	sn_date = []
	snj_date=0
	for i in range(len(file_names)):
		date=""
		name=""
		for j in range(len(sn_name)):		
			if sn_name[j] in file_names[i]:
				date=file_names[i]
				name=sn_name[j]
				if jdate[j] != 99999.9:
					snj_date = jdate[j]
		x = 24 + 2*len(name)	#18(directory length) +1 (/) + 2 (sn) + len(name)-1 +2 (//) + 2 (sn) + len(name) -1 + 1 (-)
		date = date[x:]
		date = date[:8] #length of date string
		sn_date.append(date)

	max_delta=10000000
	max_light=''
	for i in range(len(sn_date)):
		date=sn_date[i]
		year=date[:4]
		month=date[4:-2]
		day=date[6:]
		julian_date=float(year)*365+(float(month)-1)*365/12+float(day)
		if (julian_date-snj_date)<max_delta:
			max_delta = julian_date-snj_date
			max_light = file_names[i]
	
	return max_light

spectra_files=[]
for i in range(len(sn_folders)):
	spectra_files.append(find_max(glob.glob(sn_folders[i]+"/*.flm")))

ascii.write([spectra_files], 'MaxSpectra.flm', names=['col1'])