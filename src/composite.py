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
from scipy.optimize import curve_fit

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
#change this address to whereever you locally stored the SNe.db
con = sq3.connect('../../../SNe.db')
cur = con.cursor()
def grab(sql_input, Full_query):
    SN_Array = []
    cur.execute(sql_input)
    #at some point this should be more modular but for now I'm only going to accept the full query
    for row in cur:
        if sql_input == Full_query:
            SN = supernova()
            SN.filename = row[0]
            SN.name = row[1]
            SN.redshift = row[2]
	    SN.phase = row[3]
            SN.minwave = row[4]
            SN.maxwave = row[5]
	    SN.SNR = row[10]
            #spectra = msg.unpackb(row[7])
            #SN.spectrum = spectra
	    interp = msg.unpackb(row[12])
	    SN.interp = interp
	    try:
		SN.wavelength = SN.interp[0,:]
		SN.flux = SN.interp[1,:]
		SN.variance = SN.interp[2,:]
	    except TypeError:
		continue
	    full_array.append(SN)
            SN_Array.append(SN)
            #print SN.interp
	else:
	    print "Invalid query...more support will come"
    print len(SN_Array)
    #cut the array down to be more manageable
    SN_Array = SN_Array[0:50]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'variance')]
    return SN_Array
    
"""
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
"""
#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def makearray(SN_Array):
	fluxes = []
	errors = []
	flux = []
	error = []
	for SN in SN_Array:
	    #doesn't need to be truncated if data is interpolated and aligned
	    flux = SN.flux
	    error = SN.variance
	    wavelength = SN.wavelength
	    red = SN.redshift
	    #age = SN.ages[lowindex:highindex]
	    if len(fluxes) == 0:
			fluxes = np.array([flux])
			errors = np.array([error])
			reds = np.array([red])
			#ages = np.array([age])
	    else:
			try:
				fluxes = np.append(fluxes, np.array([flux]),axis=1)
				errors = np.append(errors, np.array([error]), axis=1)
				reds = np.append(reds, np.array([red]), axis = 1)
				#ages = np.append(ages, np.array([age]), axis = 0)
			except ValueError:
				continue
	return fluxes, errors

    


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
    
#Here's the function that scales spectra based on the most recent composite. It gets run multiple times if there are non-overlapping spectra.
def scfunc(x,a):
    return a*x

def fluxscale(tempflux, flux, error, lowindex, highindex):
    # Parameters:
    #     tempflux = template flux array
    #     flux     = spectrum array needs to be scaled
    #     error    = error of the spectrum needs to be scaed, use the inverse as the weighting.
    #     lowindex, highindex = wavelength range to be used to get the scale factor
    scale = curve_fit(scfunc,flux[lowindex:highindex],tempflux[lowindex:highindex])
    return scale
    
def scale(fluxes, errors, compare, wave, SN, lowindex, highindex):
    #scale to that one initialy, then scale to the current composite
    factors = []
    for i in xrange(len(wave)):
        try:
            if round(wave[i]) == round(SN.wavelength[i]):
                try:
                    factors.append(compare[i] / SN.flux[i])
                except IndexError:
                    continue
        except IndexError:
            continue
    
    scale_factor = fluxscale(compare, fluxes, errors, lowindex, highindex)
    #SN.flux[lowindex:highindex] *= scale_factor
    #SN.error[lowindex:highindex] *= scale_factor
    #scale_factor = np.mean(factors)
    SN.flux *= scale_factor
    SN.variance *= scale_factor**-2
    #plt.subplot(311)
    #plt.plot(SN.wavelength, SN.flux)
    print "Spectrum", SN.name, "scaled at factor", scale_factor
    return SN, scale_factor
    
    
    
#averages with weights based on the given errors in .flm files
def average(fluxes, errors):
	#redshift stuff wasn't working so we aren't dealing with it right now
	#avg_red = compare_spectrum.redshifts
	#redshifts = compare_spectrum.redshifts
	avg_flux = np.average(fluxes, weights = 1.0/errors, axis=0)
	avg_red = np.average(reds, weights = 1.0/errors, axis = 0)
	avg_age = np.average(ages, weights = 1.0/errors, axis = 0)
	compare_spectrum = SN_Array[0]
	compare_spectrum.flux = avg_flux
	compare_spectrum.redshifts = avg_red
	compare_spectrum.ages = avg_ages
	# Add residual formula?
	return compare_spectrum
	return avg_flux
def main():
    SN_Array = []
    #Accept SQL query as input and then grab what we need
    print "Query format: SELECT [items] FROM Supernovae"
    print "Optional at the end: ORDER BY [attribute] DESC"
    Full_query = "SELECT * FROM Supernovae"
    print "Full Query:", Full_query
    #sql_input = str(raw_input("Enter a SQL Query---> "))
    sql_input = Full_query
    SN_Array = grab(sql_input, Full_query)
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
    good = np.where(len(np.where(((wavemin>wmin) & (wavemax<wmax))>50))) #& (SN.SNR>.8*max(SN.SNR)))
    template=supernova()
    template.flux=np.array([composite.flux[good]])
    template.wavelength=np.array([composite.wavelength[good]])
    zeros=1
    tempzeros=0
    while (zeros!=tempzeros):
        for SN in SN_Array:
            SN, wavemin, wavemax, lowindex, highindex = cut(composite, SN, SN_Array, wavemin, wavemax)
        scales=[]
	fluxes, errors = makearray(SN_Array)
	print fluxes, errors
        for SN ,i in zip(SN_Array, xrange(len(SN_Array))):
	    SN,scale_factor = scale(fluxes[i,:], errors[i,:], fluxes[0,:], template.wavelength, SN, lowindex, highindex)
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
    c_file=str(raw_input("Create a file for data? (y/n)"))
    #if c_file=='y':
	#	f_name='composite,'+min(composite.phases)+'.'+max(composite.phases)+'.'+min(composite.redshifts)+'.'max(composite.redshifts)+'...--'+np.average(composite.phases)+'.'+np.average(composite.redshifts)+len(SN_Array)+'SN'
	#	#phase_min.phase_max.deltam15_min.deltam15_max. ... --avg_phase.avg_deltam15... --#SN
	#	table.write(f_name,format='ascii')
    #else:
	#	return table
