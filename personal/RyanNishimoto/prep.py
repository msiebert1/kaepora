import os
import glob
#import astropy
import numpy as np
import scipy.interpolate as intp

#list of files
spectra_files = glob.glob ('../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []

junk_data = []

#number of spectra to modify
num = 20

#get data, pathnames
for i in range(num):
	try:
		spectra_data.append(np.loadtxt(spectra_files[i]))
		file_path.append(spectra_files[i][14:-4])
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)

#table containing sn names, redshifts, etc.
sn_parameters = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)

#holds sn name
sn = []
#holds redshift value
z = []
#holds B-V
bv = []

#get relevent parameters needed for calculations
for i in range(len(sn_parameters)):
	sn.append(sn_parameters[i][0])
	z.append(sn_parameters[i][1])
	bv.append(sn_parameters[i][10])

print file_path
#deredden the spectra 
"""
get extinction in magnitudes at wavelengths

inverse flux transmission fraction IFFT == 10^(.4*extinction)
IFFT * flux value -->deredden spectra
"""
#deredshift the spectra
for i in range(num):#go through selected spectra data
	for j in range(len(sn)):#go through list of SN parameters
		if sn[j] in file_path[i]:#SN with parameter matches the path
			#print "\n",sn[j]			
			#print spectra_data[i][:,0]
			#print "redshift",z[j]
			spectra_data[i][:,0] /= (1+z[j])			
			#print spectra_data[i][:,0]



