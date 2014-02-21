import numpy as np
import os
from astropy.table import Table
from astropy.io import ascii
import scipy.interpolate as inter
import glob
import matplotlib.pyplot as plt
import math

spectra_folders = glob.glob("../../../data/cfa/*/*.flm")
sn_files=[]
bad_files=[]

#reads the SN folders and excludes the other files
for i in range(len(spectra_folders)):	
	if '.dat' in spectra_folders[i]:
		bad_files.append(spectra_folders[i])
	elif 'snf' in spectra_folders[i]:
		bad_files.append(spectra_folders[i])
	else:
		sn_files.append(spectra_folders[i])
		
parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
mjd = Table.read('../../../data/cfa/cfasnIa_mjdspec.dat', format='ascii')
sn_file = mjd["col1"]
sn_mjd=mjd["col2"]
sn_name = parameters["col1"]
jdate=parameters["col3"]
delta=[]

for i in range(len(sn_files)):
	file_name=sn_files[i]
	date=0
	max_date=0
	for j in range(len(sn_file)):
		if sn_file[j] in file_name:
			date=sn_mjd[j]
	for k in range(len(sn_name)):
		if sn_name[k] in file_name:
			max_date=jdate[k]
	delta.append(date-max_date)

#writes delta t spectrum files to a flm file
max_spectra=Table([sn_files,delta],names=('File','delta_t'))
max_spectra.write('SpectraTime.dat',format='ascii')
