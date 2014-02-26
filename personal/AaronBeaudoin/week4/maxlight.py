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

#reads the SN folders and excludes the other files
for i in range(len(spectra_folders)):
	if '.dat' in spectra_folders[i]:
		bad_folders.append(spectra_folders[i])
	elif 'snf' in spectra_folders[i]:
		bad_folders.append(spectra_folders[i])
	else:
		sn_folders.append(spectra_folders[i])

#Finds the max light spectrum of one SN by comparing the date the spectrum was taken and the MJD given in the parameter file
def find_max(file_names):
	parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
	mjd = Table.read('../../../data/cfa/cfasnIa_mjdspec.dat', format='ascii')
	sn_file = mjd["col1"]
	sn_mjd=mjd["col2"]
	sn_name = parameters["col1"]
	jdate=parameters["col3"]
	sn_date = []
	spectra_name=[]
	snj_date=0
	#finds MJD for the SN by searching the parameter file
	for i in range(len(file_names)):
		date=""
		name=""
		for j in range(len(sn_name)):		
			if sn_name[j] in file_names[i]:
				if jdate[j] != 99999.9:
					snj_date = jdate[j]
		for k in range(len(sn_file)):
			if sn_file[k] in file_names[i]:
				date = sn_mjd[k]
		sn_date.append(date)
	
	#compares the MJD with the date of each spectra, and returns the file name of the correct one
	max_delta=10000000
	max_light=''
	for i in range(len(sn_date)):
		date=sn_date[i]
		if (date-snj_date)<max_delta:
			max_delta = date-snj_date
			max_light = file_names[i]
	
	return max_light
	
def find_name(spectra_files):
	parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
	sn_name = parameters["col1"]
	spectra_names=[]
	for i in range(len(spectra_files)):
		for j in range(len(sn_name)):
			if sn_name[j] in spectra_files[i]:
				spectra_names.append(sn_name[j])
	return spectra_names

#finds the max light spectrum for every SN

spectra_files=[]
for i in range(len(sn_folders)):
	spectra_files.append(find_max(glob.glob(sn_folders[i]+"/*.flm")))
	
spectra_names=find_name(spectra_files)

parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
mjd = Table.read('../../../data/cfa/cfasnIa_mjdspec.dat', format='ascii')
sn_file = mjd["col1"]
sn_mjd=mjd["col2"]
sn_name = parameters["col1"]
jdate=parameters["col3"]
delta=[]

for i in range(len(spectra_files)):
	file_name=spectra_files[i]
	date=0
	max_date=0
	for j in range(len(sn_file)):
		if sn_file[j] in file_name:
			date=sn_mjd[j]
	for k in range(len(sn_name)):
		if sn_name[k] in file_name:
			max_date=jdate[k]
	delta.append(date-max_date)


#writes max light spectrum files to a flm file
max_spectra=Table([spectra_names,spectra_files,delta],names=('col1','col2','col3'))
max_spectra.write('MaxSpectra.dat',format='ascii')

#ascii.write([spectra_names,spectra_files], 'MaxSpectra.flm', names=['col1','col2'])

