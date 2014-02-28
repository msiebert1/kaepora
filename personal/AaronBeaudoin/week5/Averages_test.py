#!/usr/bin/env python -i
'''
Created on Feb 14, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table

class supernova(object):
    """Attributes can be added"""


#initializes things and grabs the file lists
print "Reading max light spectra from files..."
SN_Array = []
compare_spectrum = []
file_list = []
file_list = glob.glob("../../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../week4/MaxSpectra.dat", dtype = 'str', delimiter = " ", skiprows = 1)
#connect to the database, may be changed later
con = sq3.connect('../../MichaelSchubert/SNe.db')
cur = con.cursor()
names = []
j = 1
#Reads the data of max light spectra so only one per supernova is used
for line in max_light[0:50]:
	SN = supernova()
	SN.address = line[1]
	SN.name = line[0]
	data = np.loadtxt(SN.address)
	SN.wavelength = data[:,0]
	SN.flux = data[:,1]
	SN.residual = np.zeros(len(data[:,1]))   
	SN.age = line[2]
	error = data[:,2]
	if not all(x == 0.0 for x in error):
		SN.error = error
	SN_Array.append(SN)
	print j
	j += 1
print j, 'supernovae fluxed'
#Reads from database, may be changed later
print "Reading properties from database..."
j = 0
for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
	if row[0] == SN.address:
	    SN.filename = row[0]
	    SN.redshifts = row[2]
	    SN.minwave = row[3]
	    SN.maxwave = row[4]
	    names.append(SN.name)
	    print j
	    j += 1
	else:
	    continue
print len(SN_Array), "items found"
#Anything with blank error columns gets removed, since we can't make a weighted average with it.
#This can be changed to just put them in a separate array and then do a non-weighted average...if that's something we want
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]
print "Done checking. ", len(SN_Array), "items remain"       

#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

bad_range_Array = []

#averages with weights based on the given errors in .flm files
def average(SN_Array):
	#redshift stuff wasn't working so we aren't dealing with it right now
	#avg_red = compare_spectrum.redshifts
	#redshifts = compare_spectrum.redshifts
	lowindex = 0
	highindex = 2000
	fluxes = []
	errors = []
	flux = []
	error = []
	for SN in SN_Array:
	    #doesn't need to be truncated if data is interpolated and aligned
	    flux = SN.flux[0:2000]
	    error = SN.error[0:2000]
	    wavelength = SN.wavelength[0:2000]
	    if len(fluxes) == 0:
		fluxes = np.array([flux])
		errors = np.array([error])
	    else:
		fluxes = np.append(fluxes, np.array([flux]),axis=0)
		errors = np.append(errors, np.array([error]), axis=0)
	avg_flux = np.average(fluxes, weights = 1.0/errors, axis=0)
        compare_spectrum = SN_Array[0]
	compare_spectrum.flux = avg_flux
	#compare_spectrum.redshifts = avg_red
        #compare_spectrum.ages = avg_ages
	# Add residual formula?
	return compare_spectrum
#Here's the function that scales spectra based on the most recent composite. It gets run multiple times if there are non-overlapping spectra.
def splice(compare_spectrum,SN_Array,l_wave,h_wave):
	q = 1
	new_array=[]
	for SN in SN_Array:
		low_wave = l_wave
		high_wave = h_wave
	#     plt.plot(SN.wavelength, SN.flux,label=SN.name)
		if np.min(SN.wavelength)>= low_wave:
			low_wave = np.min(SN.wavelength)
		if np.max(SN.wavelength) <= high_wave:
			high_wave = np.max(SN.wavelength)
		if (np.min(SN.wavelength) >= high_wave) | (np.max(SN.wavelength) <= low_wave):
			new_array.append(SN)
		#checks to see if there is any overlap, then creates a separate array to be dealt with later
		temp = np.abs(SN.wavelength-low_wave)
		lowindex = np.where(temp == np.min(temp))
		temp = np.abs(SN.wavelength-high_wave)
		highindex = np.where(temp == np.min(temp))
		lowindex = lowindex[0]
		highindex = highindex[0]
		low = SN.wavelength[lowindex]
		high = SN.wavelength[highindex]
		SN.flux /= np.median(SN.flux)
	#    plt.plot(SN.wavelength, SN.flux,label=SN.name)
		print lowindex, "low index", highindex, "high index"
		factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
		scale_factor = np.mean(factors)
		SN.flux[lowindex:highindex] *= scale_factor
		SN.error[lowindex:highindex] *= scale_factor
		#plt.subplot(311)
		#plt.plot(SN.wavelength, SN.flux)
		print SN.address
		print "Spectrum", q, "scaled at factor", scale_factor
		q += 1
	compare_spectrum = average(SN_Array)
	return compare_spectrum, new_array

"""scale"""
print "Scaling.."

#checks for any non-overlapping spectra, then expands range if necessary
i=0
spectrum=SN_Array[0]
array=SN_Array
min=3000
max=7000
while (i==0):
	spectrum,array=splice(spectrum,array,min,max)
	if (len(array) == 0):
		i+=1
	else:
		min=min-100
		max+=100
		
compare_spectrum=spectrum

print compare_spectrum.flux, "Composite Spectrum Scaled Flux"
print compare_spectrum.error, "Composite Spectrum Error"
    	
#This makes plots and saves the composite
print len(compare_spectrum.wavelength), len(compare_spectrum.flux)
"""composite"""
composite_file = Table([compare_spectrum.wavelength[0:2000], compare_spectrum.flux[0:2000], compare_spectrum.error[0:2000]], names = ('Wavelength', 'Scaled_Flux', 'Error'))
composite_file.write('CompositeSpectrum.dat', format='ascii')
plt.figure(1)
plt.subplot(211)
plt.plot(compare_spectrum.wavelength[0:2000], compare_spectrum.flux[0:2000],label='Weighted Composite')
plt.xlabel('Wavelength (A)')
plt.ylabel('Scaled Flux')
plt.xlim(4000,7000)
#plt.ylim(0,2.8)
plt.subplot(212)
plt.plot(compare_spectrum.wavelength[0:2000], compare_spectrum.error[0:2000], label='Error')
plt.ylabel('Error')
#plt.plot(compare_spectrum.wavelength, compare_spectrum.residual, label='Residual')
plt.xlim(4000,7000)
#plt.subplot(313)
#plt.plot(compare_spectrum.wavelength, compare_spectrum.redshifts, label='Redshift')
#plt.ylabel('Avg. Redshift')
#plt.subplot(414)
#plt.plot(compare_spectrum.wavelength, compare_spectrum.ages, label='Age')
#plt.ylabel('Age')
plt.subplot(211)
plt.legend(prop={'size':10})
plt.savefig('Composite.png')
plt.show()