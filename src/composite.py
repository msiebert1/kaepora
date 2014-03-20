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
from scipy.optimize import curve_fit, minimize

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
    
class Parameters:
    """Not sure what goes here"""

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
            SN          = supernova()
            SN.filename = row[0]
            SN.name     = row[1]
            SN.redshift = row[2]
	    SN.phase    = row[3]
            SN.minwave  = row[4]
            SN.maxwave  = row[5]
	    SN.SNR      = row[10]
            #spectra = msg.unpackb(row[7])
            #SN.spectrum = spectra
	    interp      = msg.unpackb(row[12])
	    SN.interp   = interp
	    try:
		SN.wavelength = SN.interp[0,:]
		SN.flux       = SN.interp[1,:]
		SN.variance   = SN.interp[2,:]
		SN.ivar = 1/(SN.variance**2)
	    except TypeError:
		continue
	    full_array.append(SN)
            SN_Array.append(SN)
            #print SN.interp
	else:
	    print "Invalid query...more support will come"
    print len(SN_Array), "spectra found"
    #cut the array down to be more manageable
    SN_Array = SN_Array[0:100]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'variance')]
    print len(SN_Array), "spectra remain"
    for SN in SN_Array[0:5]:
	plt.plot(SN.wavelength, SN.flux)
	plt.plot(SN.wavelength, 1/(SN.ivar ** .5))
	plt.show()
	plt.savefig(SN.name + ' test spectrum.png')
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
		low2 = np.where(waves[i]==find_nearest(waves[i], round(SN_Array[i].minwave)))
		high2 = np.where(waves[i]==find_nearest(waves[i], round(SN_Array[i].maxwave)))
		print low2[0], high2[0]
		#print fluxes[i,:][low2[0]:high2[0]]
		factors = np.array(fluxes[0,:][(low2[0]+1):high2[0]]) / np.array(fluxes[i,:][(low2[0]+1):high2[0]])
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

def scale_func(params, in_data, out_data, error):
    scale = params

    model = scale * in_data
    #previously, this was returning an array, and minimize didn't like that...so I made it a single value like this.
    #I hope that's reasonable
    return np.mean((out_data - model)/error)

def find_scales(SN_Array, temp_flux, temp_ivar):
    scales = []
    print "Finding scales..."
    #loop over each SN in the array
    for SN in SN_Array:
        #grab out the flux and inverse variance for that SN
        flux = SN.flux
        ivar = SN.ivar
	tempflux = temp_flux
        #Make the combined inverse variance function.  Zeros should multiply to get zeros
        overlap = temp_ivar * ivar
        n_overlap = len([x for x in overlap if x > 0])
	#print n_overlap

        if n_overlap < 100:

            #If there is insufficient overlap, the scale is zero.
            scales = np.append(scales, np.array([0]), axis = 0)

        else:

            #Otherwise, fit things
            params = Parameters()
            params.scale = 1.0
	    #print params.scale
            #Find the appropriate values for scaling
            good = np.where(overlap > 0)
	    flux = np.array([flux[good]])
	    tempflux = np.array([temp_flux[good]])
	    ivar = np.array([ivar[good]])
            result = minimize(scale_func, params.scale, args=(flux, tempflux, ivar))
	    params.scale = float(result.x)
	    print "Scale factor = ", params.scale
            #Put the fitted value in the array
#            scales = np.append(scales, np.array([result]), axis = 0)
            scales = np.append(scales, np.array([params.scale]), axis = 0)

    return scales

def scale_data(SN_Array, scales):
    print "Scaling..."
    for i in range(len(SN_Array)):
	SN_Array[i].flux *= scales[i]
	ivar1 = SN_Array[i].ivar
	ivar2 = 1/scales[i]**2 * ivar1
	SN_Array[i].ivar = ivar2
    return SN_Array

#averages with weights based on the given errors in .flm files
def average(SN_Array, template):
	print "Averaging..."
	#print fluxes, errors
	fluxes = []
	errors = []
	for SN in SN_Array:
	    if len(fluxes) == 0:
		fluxes = np.array([SN.flux])
		errors = np.array([SN.ivar])
		#waves = np.array([wavelength])
 		#reds = np.array([red])
		#ages = np.array([age])
	    else:
		try:
		    fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
		    errors = np.append(errors, np.array([SN.ivar]), axis=0)
		    #waves = np.append(waves, np.array([wavelength]), axis=0)
		    #reds = np.append(reds, np.array([red]), axis = 0)
		    #ages = np.append(ages, np.array([age]), axis = 0)
		except ValueError:
		    print "oh god what is happening"
	for i in range(len(SN_Array)):
	    for j in range(len(fluxes[0,:])):
		if np.isnan(fluxes[i,j]):
		    fluxes[i,j] = 0
		    errors[i,j] = 0
	#print fluxes, errors
	#print errors[:,np.where(errors[0,:]!=0)]
	template.flux = np.average(fluxes[:,np.where(errors[0,:]!=0)], weights = errors[:,np.where(errors[0,:]!=0)], axis=0)
	template.flux = template.flux[0]
	template.ivar = np.average(errors[:,np.where(errors[0,:]!=0)], weights = errors[:,np.where(errors[0,:]!=0)], axis=0)
	template.ivar = template.ivar[0]
	template.wavelength = template.wavelength[np.where(errors[0,:]!=0)]
	print template.flux, template.ivar, template.wavelength
	return template

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
    print good
    template = supernova()
    template = SN_Array[good[0]]
    template = composite
    
    i = 0
    n_start = 0
    n_end = 1
    while (n_start != n_end):
	# Instead of trimming things, I think a better approach would be to make an array that is the minimum of the inverse variances of the template and comparison spectrum.  That will be zero where they don't overlap.  Then you just select the indices corresponding to non-zero values.  No for loops necessary.
        scales=[]
	n_start = len([x for x in scales if x>0])
        
	#waves, fluxes, errors, reds, factor = makearray(SN_Array)
	#scales.append(factor)
        
	scales = find_scales(SN_Array, template.flux, template.ivar)
	n_scale = len([x for x in scales if x>0])
	SN_Array = scale_data(SN_Array, scales)

	#low_overlap, SN_Array = overlap(waves, SN_Array)
	#print len(low_overlap), "do not overlap.", len(SN_Array), "spectra being averaged."
	#I think that you can scale things in the makearray function.  That would make things a little more efficient and cleaner.
	#Below can be done cleaner and without a for loop.

        template = average(SN_Array, template)
        n_end = n_scale
	n_start = n_end
	
    print "Done."
    plt.plot(template.wavelength, template.flux)
    plt.savefig('Test Composite.png')
    plt.show()
    #Either writes data to file, or returns it to user
    table = Table([template.wavelength, template.flux, template.ivar], names = ('Wavelength', 'Flux', 'Variance'))
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
