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

SN_Array = []
compare_spectrum = []

#names = np.loadtxt("", usecols = (0)) #"need to grab all the other data too"
#for row in names:
#    SN = supernova()
#    SN.input = row[0]
#    SN_Array.append(SN)
#print len(names), #"supernovae found"
file_list = []
file_list = glob.glob("../../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../week4/MaxSpectra.dat", dtype = 'str', delimiter = " ", skiprows = 1)

con = sq3.connect('../../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading max light spectra from files..."
names = []
j = 1
<<<<<<< HEAD
for line in max_light[0:290]:
	SN = supernova()
=======
for line in max_light[0:50]:
        SN = supernova()
>>>>>>> 9169c6c8ebf1d973c67a835a6b451c027039aeec
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

#SN_Array = SN_Array[0:20]
print len(SN_Array), "items found"
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]
print "Done checking. ", len(SN_Array), "items remain"       

#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

bad_range_Array = []

#averages with weights based on the given errors in .flm files
def average(compare_spectrum,SN):
	if SN == SN_Array[0]:
		compare_spectrum == SN
		return compare_spectrum
	avg_flux = compare_spectrum.flux
	#avg_red = compare_spectrum.redshifts
	mean_flux = compare_spectrum.flux
	#redshifts = compare_spectrum.redshifts
	lowindex = find_nearest(compare_spectrum.wavelength,np.min(SN.wavelength))
	highindex = find_nearest(compare_spectrum.wavelength,np.max(SN.wavelength))
#should be np.where(compare_spectrum.wavelength == np.min(SN.wavelength) if data aligned)
	for i in lowindex+np.arange(highindex-lowindex):
#        avg_flux[i] = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
		try:
			if ((SN.error[i] != 0) & (compare_spectrum.error[i] != 0)):
				fluxes = np.array([compare_spectrum.flux[i],SN.flux[i]])
				weights = np.array([1./compare_spectrum.error[i], 1./SN.error[i]])
				#redshifts = np.array([compare_spectrum.redshifts,SN.redshifts])
				#ages = np.array([compare_spectrum.age[i],SN.age[i]])
				avg_flux[i] = np.average(fluxes,weights=weights)
				#avg_red[i] = np.average(redshifts, weights = weights) 
				#compare_spectrum.error[i] = math.sqrt((compare_spectrum.flux[i]-avg_flux[i])**2 + (SN.flux[i]-avg_flux[i])**2)
				#compare_spectrum.redshifts[i] = np.average(redshifts, weights=weights)
				#avg_ages[i] = np.average(ages, weights=weights)
				compare_spectrum.error[i] = math.sqrt((weights[0]**2)*(compare_spectrum.error[i])**2 + (weights[1]**2)*(SN.error[i])**2)
			else:
				break
		except IndexError:
			print "No flux data?"
			bad_range_Array.append(SN)
			break
            
	compare_spectrum.flux = avg_flux
	#compare_spectrum.redshifts = avg_red
        #compare_spectrum.ages = avg_ages
# Add residual formula?
	return compare_spectrum

"""scale"""
print "Scaling.."

#compare_spectrum = SN_Array[0]
#scales, averages, weighted average
<<<<<<< HEAD
#plt.figure(1)
#out_of_range_low = []
#out_of_range_high = []

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
		if (min(SN.wavelength) > high_wave) | (max(SN.wavelength) < low_wave):
			new_array.append(SN)
		#check the number of scaled versus the number before scaling, if != then expand range and run again
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
		print q, "items scaled at factor", scale_factor
		compare_spectrum = average(compare_spectrum,SN)
		q += 1
	return compare_spectrum, new_array

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
	
print compare_spectrum.error
=======
plt.figure(1)
out_of_range_low = []
out_of_range_high = []
for SN in SN_Array:
    low_wave = 3000
    high_wave = 7000
#     plt.plot(SN.wavelength, SN.flux,label=SN.name)
    if np.min(SN.wavelength)>= low_wave:
        low_wave = np.min(SN.wavelength)
    if np.max(SN.wavelength) <= high_wave:
        high_wave = np.max(SN.wavelength)
    #check the number of scaled versus the number before scaling, if != then expand range and run again
    temp = np.abs(SN.wavelength-low_wave)
    lowindex = np.where(temp == np.min(temp))
    temp = np.abs(SN.wavelength-high_wave)
    highindex = np.where(temp == np.min(temp))
    lowindex = lowindex[0]
    highindex = highindex[0]
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
    print lowindex, "low index", highindex, "high index"
    factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
    scale_factor = np.mean(factors)
    SN.flux[lowindex:highindex] *= scale_factor
    SN.error[lowindex:highindex] *= scale_factor
    #plt.subplot(311)
    #plt.plot(SN.wavelength, SN.flux)
    print "Spectrum", q, "scaled at factor", scale_factor
    compare_spectrum = average(compare_spectrum,SN)
    q += 1
print compare_spectrum.flux, "Composite Spectrum Scaled Flux"
print compare_spectrum.error, "Composite Spectrum Error"
>>>>>>> 9169c6c8ebf1d973c67a835a6b451c027039aeec
    	
        
print len(compare_spectrum.wavelength)
"""composite"""
composite_file = Table([compare_spectrum.wavelength, compare_spectrum.flux, compare_spectrum.error], names = ('Wavelength', 'Scaled_Flux', 'Error'))
composite_file.write('CompositeSpectrum.dat', format='ascii')
plt.figure(1)
plt.subplot(211)
plt.plot(compare_spectrum.wavelength, compare_spectrum.flux,label='Weighted Composite')
plt.xlabel('Wavelength (A)')
plt.ylabel('Scaled Flux')
plt.xlim(4000,7000)
#plt.ylim(0,2.8)
plt.subplot(212)
plt.plot(compare_spectrum.wavelength, compare_spectrum.error, label='Error')
plt.ylabel('Error')
plt.plot(compare_spectrum.wavelength, compare_spectrum.residual, label='Residual')
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