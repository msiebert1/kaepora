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
from scipy.optimize import leastsq
import file_name

np.set_printoptions(threshold=np.nan)
mn.patch()

#Sets up some lists for later
SN_Array = []
full_array = []
compare_spectrum = []
#This only gets used if selecting by max light spectra
#max_light = []
#max_light = np.loadtxt("../personal/AaronBeaudoin/week4/MaxSpectra.dat", dtype = 'str', delimiter = " ", skiprows = 1)

class supernova(object):
    """Attributes can be added"""
    
class Parameters:
    """Not sure what goes here"""

#Connect to database
#change this address to whereever you locally stored the SNe.db
con = sq3.connect('../../../SNe.db')
#con = sq3.connect('../../temp/SNe.db')
cur = con.cursor()

#Pulls in all columns from the database for the selected query
def grab(sql_input, Full_query):
    print "Collecting data..."
    SN_Array = []
    cur.execute(sql_input)
    for row in cur:
        if sql_input == Full_query:
            SN           = supernova()
            SN.filename  = row[0]
            SN.name      = row[1]
	    SN.source    = row[2]
            SN.redshift  = row[3]
	    SN.phase     = row[4]
            SN.minwave   = row[5]
            SN.maxwave   = row[6]
	    SN.dm15      = row[7]
	    SN.m_b       = row[8]
	    SN.B_minus_v = row[9]
	    SN.targeted  = row[10]
	    SN.SNR       = row[11]
	    interp       = msg.unpackb(row[12])
	    SN.interp    = interp
	    try:
		SN.wavelength = SN.interp[0,:]
		SN.flux       = SN.interp[1,:]
		SN.ivar       = SN.interp[2,:]
	    except TypeError:
		continue
	    full_array.append(SN)
            SN_Array.append(SN)
	else:
	    print "Invalid query...more support will come"
    print len(SN_Array), "spectra found"

    #cut the array down to be more manageable
    #Used mostly in testing, if you want the full array of whatever you're looking at, comment this line out
    #SN_Array = SN_Array[0:100]
    
    for SN in SN_Array:
	for i in range(len(SN.flux)):
	    if np.isnan(SN.flux[i]):
		SN.flux[i] = 0
	    if np.isnan(SN.ivar[i]):
		SN.ivar[i] = 0
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'ivar')]
    print len(SN_Array), "spectra remain"
    return SN_Array


"""
#This is something that could be implemented in the future
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

#This is the model for how scales should be applied, used in the find_scales function
def scale_func(vars, in_data, out_data):
    
    scale = vars[0]
    model = scale * in_data
    output = model
    return output[:,0]

#...Finds scales
def find_scales(SN_Array, temp_flux, temp_ivar):
    min_overlap = 300
    scales = []
    print "Finding scales..."
    #loop over each SN in the array
    for SN in SN_Array:
        #grab out the flux and inverse variance for that SN
        flux = SN.flux
        ivar = SN.ivar
        overlap = temp_ivar * ivar
        n_overlap = len([x for x in overlap if x > 0])
	
        if n_overlap < min_overlap:

            #If there is insufficient overlap, the scale is zero.
            scales = np.append(scales, np.array([0]), axis = 0)

        else:
            #Otherwise, fit things
            vars = [1.0]
            #Find the appropriate values for scaling
            good      = np.where(overlap > 0)
	    flux2     = np.array([flux[good]])
	    ivar2     = np.array([ivar[good]])
	    tempflux2 = np.array([temp_flux[good]])
            tempivar2 = np.array([temp_ivar[good]])
            totivar   = 1/(1/ivar2 + 1/tempivar2)

	    result = np.median(tempflux2/flux2)

            if result < 0:
                result = 0

	    print "Scale factor = ", result

            scales = np.append(scales, np.array([float(result)]), axis = 0)

    return scales

#Scales the data using the factors found before
badfiles = []
def scale_data(SN_Array, scales):
    print "Scaling..."
    for i in range(len(scales)):
	SN_Array[i].flux *= np.abs(scales[i])
	SN_Array[i].ivar /= (scales[i])**2
	print "Scaled at factor ", scales[i]
    return SN_Array

#averages with weights based on the given errors in .flm files
def average(SN_Array, template):
	print "Averaging..."
	#print fluxes, errors
	fluxes = []
	ivars  = []
	reds = []
	for SN in SN_Array:
	    if len(fluxes) == 0:
		fluxes = np.array([SN.flux])
		ivars  = np.array([SN.ivar])
 		reds = np.array([SN.redshift])
		phases = np.array([SN.phase])
	    else:
		try:
		    fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
		    ivars  = np.append(ivars, np.array([SN.ivar]), axis=0)
		    reds = np.append(reds, np.array([SN.redshift]), axis = 0)
		    phases = np.append(phases, np.array([SN.phase]), axis = 0)
		except ValueError:
		    print "This should never happen!"

        #Make flux/ivar mask so we can average for pixels where everything has 0 ivar
	flux_mask = np.zeros(len(fluxes[0,:]))
	ivar_mask = np.zeros(len(fluxes[0,:]))
	have_data = np.where(np.sum(ivars, axis = 0)>0)
	no_data   = np.where(np.sum(ivars, axis = 0)==0)
	ivar_mask[no_data] = 1

        #Add in flux/ivar mask
        fluxes = np.append(fluxes, np.array(flux_mask), axis=0)
        ivars  = np.append(ivars, np.array(ivar_mask), axis=0)

	has_reds  = np.where(reds != None)
	has_phase = np.where(phases != None)
	reds      = reds[has_reds]
	phases    = phases[has_phase]

	for i in range(len(fluxes)):
#	    ivars[i,:] += ivar_mask
        template.flux = np.average(fluxes, weights=ivars, axis=0)
        template.ivar = 1/np.sum(ivars, axis=0)
	template.redshift = sum(reds)/len(reds)
	#template.phase = sum(phases)/len(phases)
	#phase data isn't in the database yet
	template.ivar[no_data] = 0
	template.name = "Composite Spectrum"
	return template

def main(Full_query):
    SN_Array = []
    #Accept SQL query as input and then grab what we need
    #Full_query = "SELECT * FROM Supernovae WHERE snr > 8"
    print "SQL Query:", Full_query
    sql_input = Full_query

    SN_Array = grab(sql_input, Full_query)

    #finds the longest SN we have for our initial template
    lengths = []
    for SN in SN_Array:
        lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
    temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
    composite = temp[0]

    #scales data, makes a composite, and splices in non-overlapping data
    #Here is where we set our wavelength range for the final plot
    wmin = 4000
    wmax = 7500
    wavemin = composite.minwave
    wavemax = composite.maxwave

    #finds range of useable data
    good = np.where(len(np.where(((wavemin <= wmin) & (wavemax >= wmax)) > 100)))
    template = supernova()
    template = SN_Array[good[0]]
    template = composite
    
    #Starts our main loop
    i = 0
    n_start = 0
    n_end = 1
    scales=[]
    while (n_start != n_end):
	n_start = len([x for x in scales if x>0])
        scales=[]       
	scales = find_scales(SN_Array, template.flux, template.ivar)
	n_scale = len([x for x in scales if x>0])
	SN_Array = scale_data(SN_Array, scales)
        template = average(SN_Array, template)
        n_end = n_scale
	n_start = n_end
	
	
    print "Done."
    print "Average redshift =", template.redshift
    #print "Average phase =", template.phase
    #This next line creates a unique filename for each run based on the sample set used
    f_name = "../plots/" + file_name.make_name(SN_Array)
    template.savedname = f_name
    lowindex = np.where(template.wavelength == find_nearest(template.wavelength, wmin))
    highindex = np.where(template.wavelength == find_nearest(template.wavelength, wmax))
    plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.flux[lowindex[0]:highindex[0]])
    plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.ivar[lowindex[0]:highindex[0]])
    plt.savefig('../plots/' + f_name + '.png')
    plt.show()
    #Either writes data to file, or returns it to user
    #This part is still in progress
    table = Table([template.wavelength, template.flux, template.ivar], names = ('Wavelength', 'Flux', 'Variance'))
    c_file = str(raw_input("Create a file for data? (y/n)"))
    if c_file=='y':
		#f_name = "../plots/TestComposite"
		table.write(f_name,format='ascii.no_header')
		return template
    else:
		return template

if __name__ == "__main__":
    main()
