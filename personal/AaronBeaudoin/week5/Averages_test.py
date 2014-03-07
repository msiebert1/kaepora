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
for line in max_light[0:100]:
	SN = supernova()
	SN.address = line[1]
	SN.name = line[0]
	data = np.loadtxt(SN.address)
	SN.wavelength = data[:,0]
	SN.flux = data[:,1]
	SN.residual = np.zeros(len(data[:,1])) 
        SN.redshifts = np.zeros(len(data[:,1]))
        SN.redshifts.fill(SN.redshift)  
	SN.age = line[2]
        SN.ages = np.zeros(len(data[:,1]))
        SN.ages.fill(SN.age)
	error = data[:,2]
	if not all(x == 0.0 for x in error):
		SN.error = error
	SN_Array.append(SN)
	#print j
	j += 1
print j, 'supernovae fluxed'
#Reads from database, may be changed later
print "Reading properties from database..."
j = 0
for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
	if row[0] in SN.address:
	    SN.filename = row[0]
	    SN.redshift = row[2]
	    SN.minwave = row[3]
	    SN.maxwave = row[4]
	    names.append(SN.name)
	    #print j
	    print SN.redshift
	    j += 1
	else:
	    continue
print len(SN_Array), "items found"

"""open the velocity file, make it useable, add the needed data to the array"""
print "Reading velocities..."
v_data = np.loadtxt('../../SamRubin/foley_master_data', dtype={'names': ('name', 'redshift', 'v', 'dv'), 'formats': ('S8', 'f8', 'f8', 'f8')}, skiprows = 1, usecols = (0,1,2,3))
for SN in SN_Array:
    for row in v_data:
	if row[0] == SN.name:
	    if row[2] != -99.0:
		SN.v_si = row[2]
		SN.dv_si = row[3]
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'v_si')]
print len(SN_Array), "velocities found"

#Anything with blank error columns gets removed, since we can't make a weighted average with it.
#This can be changed to just put them in a separate array and then do a non-weighted average...if that's something we want
print "Checking for blank error columns..."
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]
print "Done checking. ", len(SN_Array), "items remain"       

#Here is some new code that might not work at all.
#The goal is to create a system which will run different comparisons between spectra based on user input
#Current separation conditions are V_SiII, host galaxy, redshift
#Takes input, runs through the python equivalent of a switch
#runs a different separation subroutine and returns two arrays, each containing a bunch of spectra that meet the criteria
#Then those two arrays get run through the compositer individually.
array1 = []
array2 = []

def split_by_v(SN_Array):
    high_v = raw_input('Velocity Boundary = ')
    high_v = float(high_v)
    for SN in SN_Array:
	if SN.v_si > high_v:
	    array1.append(SN)
	else:
	    array2.append(SN)
    return array1, array2

def split_by_host(SN_Array):
    host_type = np.float(raw_input('Input a host type --->'))
    for SN in SN_Array:
	if 1 == 1: #dummy statement until i figure out what to put there
	    array1.append(SN)
	else:
	    array2.append(SN)
    return array1, array2

def split_by_red(SN_Array):
    high_red = raw_input('Redshift Boundary = ')
    high_red = float(high_red)
    for SN in SN_Array:
	if SN.redshift > high_red:
	    array1.append(SN)
	else:
	    array2.append(SN)
    return array1, array2

def no_split(SN_Array):
    array1 = SN_Array
    array2 = []
    return array1, array2

#here's the switch
select = {"1" : split_by_v,
	  "2" : split_by_host, #doesn't work yet
	  "3" : split_by_red,
	  "4" : no_split, #doesn't do anything, just takes the input and spits it right back out
	  }
print "Do you want to split the data for comparison?"
print "1. Split by Silicon Line Velocity"
print "2. Split by Host Galaxy Type"
print "3. Split by Redshift"
print "4. No Split"
choice = raw_input('Choose a split profle --> ')
array1, array2 = select[choice](SN_Array)
print array1[0].flux
print array2[0].flux

#after this we go back into the normal composite stuff
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
            red = SN.redshifts[0:2000]
            age = SN.ages[0:2000]
	    if len(fluxes) == 0:
		fluxes = np.array([flux])
		errors = np.array([error])
                reds = np.array([red])
                ages = np.array([age])
	    else:
		fluxes = np.append(fluxes, np.array([flux]),axis=0)
		errors = np.append(errors, np.array([error]), axis=0)
                reds =  np.append(reds, np.array([red]), axis=0)
                ages =  np.append(ages, np.array([age]), axis=0)
	avg_flux = np.average(fluxes, weights = 1.0/errors, axis=0)
        avg_red = np.average(reds, weights = 1.0/errors, axis=0)
        avg_age = np.average(ages, weights = 1.0/errors, axis=0)
        compare_spectrum = SN_Array[0]
	compare_spectrum.flux = avg_flux
	compare_spectrum.redshifts = avg_red
        compare_spectrum.ages = avg_age
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
		factors = compare_spectrum.flux[0:2000] / SN.flux[0:2000]
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
min=3000
max=7000
while (i==0):
	compare_spectrum,array=splice(spectrum,SN_Array,min,max)
	spectrum1_new,array1_new=splice(spectrum,array1,min,max)
	spectrum2_new,array2_new=splice(spectrum,array2,min,max)
	
	if (len(array) == 0):
		i+=1
	else:
		min=min-100
		max+=100
		
composite1 = spectrum1_new
composite2 = spectrum2_new

#print compare_spectrum.flux, "Composite Spectrum Scaled Flux"
#print compare_spectrum.error, "Composite Spectrum Error"
    	
#This makes plots and saves the composite
print len(composite1.wavelength), len(composite2.flux)
"""composite"""
composite_file = Table([compare_spectrum.wavelength[0:2000], compare_spectrum.flux[0:2000], compare_spectrum.error[0:2000]], names = ('Wavelength', 'Scaled_Flux', 'Error'))
composite_file.write('CompositeSpectrum.dat', format='ascii')
composite_file1 = Table([composite1.wavelength[0:2000], composite1.flux[0:2000], composite1.error[0:2000]], names = ('Wavelength', 'Scaled_Flux', 'Error'))
composite_file1.write('CompositeSpectrum1.dat', format='ascii')
composite_file2 = Table([composite2.wavelength[0:2000], composite2.flux[0:2000], composite2.error[0:2000]], names = ('Wavelength', 'Scaled_Flux', 'Error'))
composite_file2.write('CompositeSpectrum2.dat', format='ascii')
plt.figure(1)
plt.subplot(211)
plt.plot(compare_spectrum.wavelength[0:2000], compare_spectrum.flux[0:2000], label='Weighted Composite')
plt.plot(composite1.wavelength[0:2000], composite1.flux[0:2000],label='Weighted Composite 1')
plt.plot(composite2.wavelength[0:2000], composite2.flux[0:2000],label='Weighted Composite 2')
plt.xlabel('Wavelength (A)')
plt.ylabel('Scaled Flux')
plt.xlim(4000,7000)
#plt.ylim(0,2.8)
plt.subplot(212)
plt.plot(composite1.wavelength[0:2000], composite1.error[0:2000], label='Error 1')
plt.plot(composite2.wavelength[0:2000], composite2.error[0:2000], label='Error 2')
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
