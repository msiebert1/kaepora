import os
import glob
from specutils import extinction as ex
from astropy.table import Table
from astropy.io import ascii
import astroquery
from astroquery.irsa_dust import IrsaDust
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter


#list of files
spectra_files = glob.glob ('../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []

junk_data = []

#number of spectra to modify
num = 1

#get data, pathnames
for i in range(num):
	try:
        	spectra_data.append(np.loadtxt(spectra_files[i]))
        	file_path.append(spectra_files[i][14:-4])
		#print file_path
             
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

#get relevent parameters needed for calculations
for i in range(len(sn_parameters)):
	sn.append(sn_parameters[i][0])
	z.append(sn_parameters[i][1])

"""
NOTE:
Use IRSA
"""
#deredden and deredshift the spectra
for i in range(num):#go through selected spectra data
	for j in range(len(sn)):#go through list of SN parameters
		if sn[j] in file_path[i]:#SN with parameter matches the path
			print "looking at SN",sn[j]
			ext = IrsaDust.get_extinction_table("SN%s"%sn[j])
			print ext[1][0],ext[1][3]
			print ext[2][0],ext[2][3]
			b = ext[1][3]
			v = ext[2][3]
						
			print "starting flux:\n",spectra_data[i][:,1]
			print "b-v value:",b,v
			spectra_data[i][:,1] = spectra_data[i][:,1]*ex.reddening(spectra_data[i][:,0],ebv=b-v,r_v=3.1,model='fm07')
			print "de-reddened flux:\n",spectra_data[i][:,1]
			print "starting wavelength:\n",spectra_data[i][:,0]
			spectra_data[i][:,0] /= (1+z[j])
			print "z:",z[j]
			print "de-red-shifted wavelength:\n",spectra_data[i][:,0]	


