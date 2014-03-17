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

np.set_printoptions(threshold=np.nan)
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
#con = sq3.connect('../../temp/SNe.db')
cur = con.cursor()
def grab(sql_input, Full_query):
    print "Collecting data..."
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
    print len(SN_Array), "spectra found"
    #cut the array down to be more manageable
    SN_Array = SN_Array[0:50]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'variance')]
    print len(SN_Array), "spectra remain"
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

def scfunc(x,a):
    return a*x

def makearray(SN_Array):
	print "Creating arrays..."
	fluxes = []
	errors = []
	flux   = []
	error  = []
	waves = []
	red = []
	reds = []

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
			waves = np.array([wavelength])
 			reds = np.array([red])
			#ages = np.array([age])
	    else:
			try:
				fluxes = np.append(fluxes, np.array([flux]), axis=0)
				errors = np.append(errors, np.array([error]), axis=0)
				waves = np.append(waves, np.array([wavelength]), axis=0)
				reds = np.append(reds, np.array([red]), axis = 0)
				#ages = np.append(ages, np.array([age]), axis = 0)
			except ValueError:
				continue
	print "Scaling..."
	for i in range(len(SN_Array)):
	    #try:
		#I still can't make curve_fit work, so I'm trying something else. -Sam
		#factor, pcov = curve_fit(scfunc, fluxes[0,:], fluxes[i,:], sigma = 1/errors[i,:])
		flux1 = np.array([value for value in fluxes[0,:] if not math.isnan(value)])
		flux2 = np.array([value for value in fluxes[i,:] if not math.isnan(value)])
		#Ideally, this only scales the region overlapping spectrum 0.
		#Somehow, using the same slice in two different arrays gives different sizes.
		#Need to fix it.
		low2 = np.where(waves[i]==find_nearest(waves[i], round(SN_Array[i].minwave)))
		high2 = np.where(waves[i]==find_nearest(waves[i], round(SN_Array[i].maxwave)))
		print low2[0], high2[0]
		factors = flux1[low2[0]:high2[0]] / flux2[low2[0]:high2[0]]
		#print factors
		factor = np.mean(factors)
		fluxes[i,:] *= factor
		errors[i,:] *= factor**-2
		print "Scaled at factor", factor
	    #except RuntimeError:
		#print "Curve-fit failed"
	#print fluxes, errors
	return waves, fluxes, errors, reds, factor

low_overlap = []
def overlap(waves, SN_Array):
    #select overlapping spectra (overlap by more than 100 angstroms)
    print "Checking for overlap..."
    common = []
    for SN in SN_Array:
        common = [val for val in SN.wavelength if val in waves[0,:]]
        if len(common) <= 100:
            low_overlap.append(SN)
            SN_Array.remove(SN)
    return low_overlap, SN_Array

#averages with weights based on the given errors in .flm files
def average(SN_Array, fluxes, errors, reds):
	print "Averaging..."
	newflux = []
	newerror = []
	for i in range(len(SN_Array)):
	    newflux[i] = np.array([value for value in fluxes[i,:] if not math.isnan(value)])
	    newerror[i] = np.array([value for value in errors[i,:] if not math.isnan(value)])
	avg_flux = np.average(newflux, weights = 1.0/newerror, axis=0)
	print avg_flux
	#avg_red = np.average(reds, weights = 1.0/errors, axis = 0)
	#avg_age = np.average(ages, weights = 1.0/errors, axis = 0)
	compare_spectrum = SN_Array[0]
	compare_spectrum.flux = avg_flux
	#compare_spectrum.redshifts = avg_red
	#compare_spectrum.ages = avg_ages
	return compare_spectrum
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
        lengths.append(len(SN.wavelength[np.where(SN.flux != 0)]))
    temp = [SN for SN in SN_Array if len(SN.wavelength) == max(lengths)]
    composite = temp[0]
#     print composite.flux

    #scales data, makes a composite, and splices in non-overlapping data
    wmin = 4500
    wmax = 6500
    wavemin = composite.minwave
    wavemax = composite.maxwave

    good = np.where(len(np.where(((wavemin > wmin) & (wavemax < wmax)) > 100))) #& (SN.SNR>.8*max(SN.SNR)))
    template = supernova()
    template.wavelength = np.array([composite.wavelength[good]])
    template.flux       = np.array([composite.flux[good]])

    i = 0
    while (i == 0):
	# Instead of trimming things, I think a better approach would be to make an array that is the minimum of the inverse variances of the template and comparison spectrum.  That will be zero where they don't overlap.  Then you just select the indices corresponding to non-zero values.  No for loops necessary.
        scales=[]
	waves, fluxes, errors, reds, factor = makearray(SN_Array)
	scales.append(factor)
	low_overlap, SN_Array = overlap(waves, SN_Array)
	print len(low_overlap), "do not overlap.", len(SN_Array), "spectra being averaged."
	#I think that you can scale things in the makearray function.  That would make things a little more efficient and cleaner.
	#Below can be done cleaner and without a for loop.
        template = average(SN_Array, fluxes, errors, reds)
        if not np.any(scales) == 0:
	    i += 1
    print "Done."
    #Either writes data to file, or returns it to user
    table = Table([template.wavelength, template.flux, template.variance], names = ('Wavelength', 'Flux', 'Variance'))
    c_file = str(raw_input("Create a file for data? (y/n)"))
    if c_file=='y':
		#f_name='composite,'+min(composite.phases)+'.'+max(composite.phases)+'.'+min(composite.redshifts)+'.'+max(composite.redshifts)+'...--'+np.average(composite.phases)+'.'+np.average(composite.redshifts)+len(SN_Array)+'SN'
		#phase_min.phase_max.deltam15_min.deltam15_max. ... --avg_phase.avg_deltam15... --#SN
		f_name = "Test Composite"
		table.write(f_name,format='ascii')
    else:
		return table

if __name__ == "__main__":
    main()
