"""
misc. notes:
read in all spectra
pixels

tweaks for composite
de-reddening
	#dust in milky way
		#using fitzpatrick law
		#find milky way reddening for particular SN (look on NED)

	#de-redshifting

	#by host galaxy (assume zero for now)
		#color of SN matters

interpolate
NO TRUNCATING
1000-20000
set inverse variance (square of error)(weight) to zero outside of wavelength ranges 

List of Files containing Spectra Data
"""
import os
import numpy as np
import astropy
from astropy.table import Table
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import glob
import specutils

#put files in an array
spectra_files = glob.glob ('../../data/cfa/*/*.flm')

#holds the wavelength/flux/weight of each spectra
spectra_data = []
#holds the pathnames of files
file_path = []
#accounts for bad files
junk_data = []

#look at only 10 files for now since my computer is super slow
for i in range(10):
	try:
		spectra_data.append(np.loadtxt(spectra_files[i]))
		file_path.append(spectra_files[14:-4])
	except ValueError:
		junk_data.append(spectra_files)
	
#table containing sn names, redshifts, etc.
sn_parameters = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
sn = []
z = []
for i in range(len(sn_parameters)):
	sn.append(sn_parameters[i][0])
	z.append(sn_parameters[i][1])

#for i in range(len(sn_parameters)):
#	print sn[i],"has redshift z=",z[i]

print specutils.extinction_fm07(spectra_data[0][0])
#for i in range(len(spectra_data)):
#	for j in range(len(sn)):
#		if sn[j] in file_path[i]:
#			spectra_data[i][0] = specutils.extinction_fm07(spectra_data[i][0])		
#			spectra_data[i][0] /= 1+z[j]
#			spectra_data[i][0] = specutils.extinction_fm07(spectra_data[i][0], 0)

