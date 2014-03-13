"""
Spectra composite program
Authors: Sam, Yixian, Aaron
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table
import msgpack as msg
import msgpack_numpy as mn
import lmfit

mn.patch()

#Sets up some lists for later
SN_Array = []
full_array = []
compare_spectrum = []
file_list = []
file_list = glob.glob("../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../personal/AaronBeaudoin/week4/MaxSpectra.dat", dtype = 'str', delimiter = " ", skiprows = 1)

class supernova(object):
    """Attributes can be added"""

#Connect to database
con = sq3.connect('../../../SNe.db')
cur = con.cursor()

#Accept SQL query as input and then grab what we need
print "Query format: SELECT [items] FROM Supernovae"
print "Optional at the end: ORDER BY [attribute] DESC"
Full_query = "SELECT Filename, SN, Redshift, MinWave, MaxWave, Spectra FROM Supernovae"
print "Full Query:", Full_query
sql_input = str(raw_input("Enter a SQL Query---> "))
cur.execute(sql_input)
#at some point this should be more modular but for now I'm only going to accept the full query
for row in cur:
    if sql_input == Full_query:
        SN = supernova()
        SN.filename = row[0]
        SN.name = row[1]
        SN.redshift = row[2]
        SN.minwave = row[3]
        SN.maxwave = row[4]
        spectra = msg.unpackb(row[5])
        SN.spectrum = spectra
        full_array.append(SN)
        SN.wavelength = SN.spectrum[:,0]
        SN.flux = SN.spectrum[:,1]
        SN.variance = SN.spectrum[:,2]
        SN.redshifts = np.zeros(len(SN.wavelength))
	SN.redshifts.fill(SN.redshift)
        #print SN.spectrum
    else:
        print "Invalid query"
print len(full_array)
#Only keeps one per supernova at max light. Condition can be changed later.
for SN in full_array:
    for row in max_light:
        if SN.filename in row[1]:
            SN_Array.append(SN)
            SN.age = row[2]
            SN.ages = np.zeros(len(SN.wavelength))
            SN.ages.fill(SN.age)
            #print SN.age
print len(SN_Array)

#cut the array down to be more manageable
SN_Array = SN_Array[0:50]
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'variance')]

"""
#create a linspace for the wavelengths
wave_min = math.floor(max(SN.minwave for SN in SN_Array))
wave_max = math.floor(min(SN.maxwave for SN in SN_Array))
waverange = np.linspace(wave_min, wave_max, (wave_max-wave_min))
waverange = np.array(waverange)
#interpolate to make my life easier
for SN in SN_Array:
    waverange = np.linspace(min(SN.wavelength), max(SN.wavelength), (max(SN.wavelength)-min(SN.wavelength)))
    SN.wavelength = np.divide(SN.wavelength, 1 + SN.redshift)
    try:
        spline1 = intp.splrep(waverange, SN.flux)
        spline2 = intp.splrep(waverange, SN.variance)
    except ValueError:
        print "Invalid data found"
    new_flux = intp.splev(waverange, spline1)
    new_flux /= np.median(new_flux)
    new_error = intp.splev(waverange, spline2)
    SN.wavelength = waverange
    SN.flux = new_flux
    SN.variance = new_error
    print SN.flux
""" 
#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def cut(compare, SN, SN_Array, min_wave, max_wave):
    #determine good wavelength range, trim spectra
    for i in xrange(len(SN.wavelength)):
        if SN.wavelength[i] == find_nearest(SN.wavelength, min_wave):
            lowindex = i
    for i in xrange(len(SN.wavelength)):
	if SN.wavelength[i] == find_nearest(SN.wavelength, max_wave):
	    highindex = i
    if highindex == 0:
        SN_Array.remove(SN)
    else:    
        print lowindex, "to", highindex
    SN.wavelength = SN.wavelength[lowindex:highindex]
    SN.flux = SN.flux[lowindex:highindex]
    SN.variance = SN.variance[lowindex:highindex]
    SN.redshifts = SN.redshifts[lowindex:highindex]
    SN.spectrum = np.array([SN.wavelength, SN.flux, SN.variance])
    return SN, min_wave, max_wave, lowindex, highindex

low_overlap = []
def overlap(compare, SN_Array):
    #select overlapping spectra (overlap by more than 100 angstroms)
    common = []
    for SN in SN_Array:
        common = [val for val in SN.wavelength if val in compare.wavelength]
        if len(common) <= 100:
            low_overlap.append(SN)
            SN_Array.remove(SN)
    return low_overlap, SN_Array
    
def scale(compare, SN,SN_Array):
    #scale to that one initialy, then scale to the current composite
    factors = []
	low_overlap, SN_Array = overlap(composite, SN_Array)
    for i in xrange(len(compare.wavelength)):
        try:
            if round(compare.wavelength[i]) == round(SN.wavelength[i]):
                try:
                    factors.append(compare.flux[i] / SN.flux[i])
                except IndexError:
                    continue
        except IndexError:
            continue
    scale_factor = np.mean(factors)
    SN.flux *= scale_factor
    SN.variance *= scale_factor**-2
    #plt.subplot(311)
    #plt.plot(SN.wavelength, SN.flux)
    print "Spectrum", SN.name, "scaled at factor", scale_factor
    return SN,scale_factor,SN_Array
    
    
    
def average(composite, SN_Array):
	fluxes = []
	errors = []
	flux = []
	error = []
	average_array = []
	waverange = []
	j=2
	for SN in SN_Array[1:]:
	    #doesn't need to be truncated if data is interpolated and aligned
	    flux = SN.flux
	    error = SN.variance
	    wavelength = SN.wavelength
	    red = SN.redshifts
	    age = SN.ages
	    for i in xrange(len(composite.wavelength)):
		composite.wavelength[i] = round(composite.wavelength[i])
	    fluxes = np.array([composite.wavelength])
	    fluxes = np.append(fluxes, np.array([composite.flux]), axis=0)
	    print fluxes
	    for i in xrange(len(wavelength)):
		wavelength[i] = round(wavelength[i])
		if composite.wavelength[i] == wavelength[i]:
		    fluxes[i] = np.append(fluxes[i], np.array([flux[i]]), axis = 1)
		    #errors = np.append(errors, np.array([error]), axis=0)
		    #reds = np.append(reds, np.array([red]), axis = 0)
		    #ages = np.append(ages, np.array([age]), axis = 0)
		else:
		    fluxes[i] = np.append(fluxes[i], np.array([0]), axis = 1)
	    j+=1
	avg_flux = np.average(fluxes, weights = 1.0/errors, axis=0)
	#avg_red = np.average(reds, weights = 1.0/errors, axis = 0)
	#avg_age = np.average(ages, weights = 1.0/errors, axis = 0)
        composite.flux = avg_flux
        
	#compare_spectrum.redshifts = avg_red
        #compare_spectrum.ages = avg_ages
	# Add residual formula?
	return compare_spectrum  
	
#finds the longest SN we have for comparison
lengths = []
for SN in SN_Array:
    lengths.append(len(SN.wavelength))
temp = [SN for SN in SN_Array if len(SN.wavelength) == max(lengths)]
composite = temp[0]
print composite.flux

#scales data, makes a composite, and splices in non-overlapping data
wmin = 4000
wmax = 6000
wavemin=composite.minwave
wavemax=composite.maxwave
good = np.where(len(np.where((wavemin>wmin) & (wavemax<wmax))>50)) #& (SN.SNR>.8*max(SN.SNR)))
template=composite.flux[good]
zeros=1
tempzeros=0

while (zeros!=tempzeros):
    for SN in SN_Array:
        SN, min_wave, max_wave, lowindex, highindex = cut(composite, SN, SN_Array, min_wave, max_wave)
    scales=[]
	for SN in SN_Array:
        SN,scale_factor,SN_Array = scale(template,SN,SN_Array)
		scales.append(scale_factor)
    template = average(template, SN_Array)
	tempzeros=0
	for i in range(len(scales)):
		if scales[i]==0:
			tempzeros+=1
	composite=template
	zeros=tempzeros
        #min_wave -= 100
        #max_wave += 100

#Either writes data to file, or returns it to user
table=Table([composite.wavelength,composite.flux,composite.variance],names=('Wavelength','Flux','Variance'))
c_file=raw_input("Create a file for data? (y/n)")
if (c_file='y'):
	f_name=raw_input("Input file name: ")
	table.write(f_name,format='ascii')
else:
	return table