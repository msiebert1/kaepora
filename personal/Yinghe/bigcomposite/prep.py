import os
import glob
from specutils import extinction as ex
import astroquery
from astroquery.ned import Ned
from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from math import floor,ceil
#import sqlite3 as sq3
#import msgpack

#con = sq3.connect('SNe.db')
#cur = con.cursor()

#list of files
spectra_files = glob.glob ('../../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []
#holds file name
file_name = []

junk_data = []

#number of spectra to modify
num = 200

#get data, pathnames
for i in range(num):
	try:
         spectra_data.append(np.loadtxt(spectra_files[i]))
         file_path.append(spectra_files[i][14:-4])
         file_name.append(spectra_files[i][27:-4])
		#print file_path
             
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)

#table containing sn names, redshifts, etc.
sn_parameters = np.genfromtxt('../../../data/cfa/cfasnIa_param.dat',dtype = None)
#table containing B and V values for determining extinction -> dereddening due to milky way
sne= np.genfromtxt('extinction.dat', dtype = None)

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
Using IRSA
"""

"""
Note: using parameters.dat file which was created from paramaters.py
parameters.py is designed to pull all relevant parameters for SNe spectra from online databases
via astroquery and place it all in a table to be pulled from later;
it would ideally do that following:
-read in all files we want to prep
-truncate file name to format "SNyear"(this is how astroquery searches for SN)
-get relevent data into table, following a format like:
SN name		Host Galaxy		Redshift	B	V	De-redden factor	CarbonPos/Neg

####
NOTE:Currently only has SN_name, B, and V values for purposes of Dereddening due to Milky way dust
####
"""

#deredden spectra to milky way
#deredshift the spectra
#deredden to host galaxy
for i in range(num):#go through selected spectra data
	for j in range(len(sn)):#go through list of SN parameters
		if sn[j] in file_path[i]:#SN with parameter matches the path
			print "\n##############################################################\nlooking at",sne[j],"matched with",sn[j]
			b = sne[j][1].astype(float)
			v = sne[j][2].astype(float)
			bv = b-v
			r = v/bv
			print "B(%s)-V(%s)=%s"%(b,v,bv)
			print "R(v) =",r
			print "starting flux:\n",spectra_data[i][:,1]
			#or use fm07 model
			#test1 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='ccm89')
			#test2 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='od94')
			spectra_data[i][:,1] *= ex.reddening(spectra_data[i][:,0],ebv = bv, r_v = 3.1, model='f99')
			print "de-reddened with specutils f99 model:\n",spectra_data[i][:,1]
			#print "de-reddened with specutils ccm89 model:\n",test1
			#print "de-reddened with specutils od94 model:\n",test2

			print "starting wavelength:\n",spectra_data[i][:,0]
			spectra_data[i][:,0] /= (1+z[j])
			print "z:",z[j]
			print "de-red-shifted wavelength:\n",spectra_data[i][:,0]
			
			#print "de-reddened by host galaxy\n",spectra_data[i][:,1]*ex.reddening(spectra_data[i][:,0],ebv = 0, r_v = r, model='f99')
			#host *= ex.reddening(spectra_data[i][:,0],ebv = bv, r_v = r, model='f99')

print "##############################################################\ndone de-reddening and de-redshifting"
##############################################################################################################################################
##############################################################################################################################################

# Data Interpolation

"""
NOTE:
This function inputs three lists containing wavelength, flux and variance.
The output will be a Table with all the fitted values.
You can change the mininum and maximum of wavelength to output, as well as pixel size in the first few lines.
For the spectra that does not cover the whole range of specified wavelength,
we output the outside values as NAN
"""


def Interpo(wave,flux,variance) :
    wave_min = 1000
    wave_max = 20000
    pix = 2
    #wavelength = np.linspace(wave_min,wave_max,(wave_max-wave_min)/pix+1)  
    wavelength = np.arange(ceil(wave_min), floor(wave_max), dtype=int, step=pix) #creates N equally spaced wavelength values
    fitted_flux = []
    fitted_var = []
    new = []
    lower = wave[0] # Find the area where interpolation is valid
    upper = wave[len(wave)-1]
    lines = np.where((wave>lower) & (wave<upper))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
    indata=inter.splrep(wave[lines],flux[lines])	#creates b-spline from new spectrum
    inerror=inter.splrep(wave[lines],variance[lines]) # doing the same with the errors
    fitted_flux=inter.splev(wavelength,indata)	#fits b-spline over wavelength range
    fitted_var=inter.splev(wavelength,inerror)   # doing the same with errors
    badlines = np.where((wavelength<lower) | (wavelength>upper))
    fitted_flux[badlines] = float('NaN')  # set the bad values to NaN !!! 
    fitted_var[badlines] = float('NaN') 
    new = Table([wavelength,fitted_flux,fitted_var],names=('col1','col2','col3')) # put the interpolated data into the new table    
    return new # return new table


# Interpolation, noise and all the stuff

from datafidelity import *  # Get variance from the datafidelity outcome

navg = [] # average noise for each spectra

for i in range(num) : 
    new_spectrum=spectra_data[i]	#declares new spectrum from list
    new_wave=new_spectrum[:,0]	#wavelengths
    new_flux=new_spectrum[:,1]	#fluxes
    var = genvar(new_wave, new_flux) #variance
#    var = new_flux*0+1
    newdata = Interpo(new_wave,new_flux,var) # Do the interpolation

    
    # Get the Noise for each spectra
    noise = new_flux/var**2.0
    navg.append(np.median(noise))

#######################################################################################################################
################### The rest is just for output testing################################################################
################### No need to implement this part to the database code ###############################################


# output data into a file 

#    output = 'testdata/modified-%s.dat'%(file_name[i])
#    ascii.write(newdata, output)
   
    # plot spectra 
    x = newdata['col1']
    y = newdata['col2']
    z = newdata['col3']
#    print z
    plt.subplot(1,2,1)
    plt.plot(x,y)
    plt.xlim(3000,7000)
    plt.subplot(1,2,2)
    plt.plot(x,z)
    plt.xlim(3000,7000)

plt.show()
plt.savefig('test_host.png')


################## Finishline of output testing #######################################################################
#######################################################################################################################


# Output of noise

#print navg
ntable = Table([file_name,navg],names=('spectra','noise'))   
ascii.write(ntable,'noise.dat')

