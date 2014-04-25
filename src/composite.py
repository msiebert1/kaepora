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
from scipy.special import erf
import file_name
###When bootstrap works right, uncomment this.
import bootstrap

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

#Connect to database
#change this address to whereever you locally stored the SNe.db
#We should all be using the same save location, since there's a .gitignore now
con = sq3.connect('../data/SNe.db')
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
            SN.velocity  = row[10]
            SN.morph     = row[11]
            SN.carbon    = row[12]
            SN.GasRich   = row[13]
            SN.SNR       = row[14]
            interp       = msg.unpackb(row[15])
            SN.interp    = interp
            try:
                SN.wavelength = SN.interp[0,:]
                SN.flux       = SN.interp[1,:]
                SN.ivar       = SN.interp[2,:]
            except TypeError:
                continue
            full_array.append(SN)
            SN_Array.append(SN)
            # print SN.filename
        else:
            print "Invalid query...more support will come"
    print len(SN_Array), "spectra found"

    #cut the array down to be more manageable
    #Used mostly in testing, if you want the full array of whatever you're looking at, comment this line out
    #SN_Array = SN_Array[0:10]
    
    #Within the interpolated spectra there are a lot of 'NaN' values
    #Now they become zeros so things work right
    for SN in SN_Array:
        SN.age = np.array(SN.flux)
        SN.dm15_array = np.array(SN.flux)
	SN.red_array  = np.array(SN.flux)
	SN.vel        = np.array(SN.flux)
        for i in range(len(SN.flux)):
            if np.isnan(SN.flux[i]):
                SN.flux[i] = 0
            if np.isnan(SN.ivar[i]):
                SN.ivar[i] = 0
            if np.isnan(SN.age[i]):
                SN.age[i]  = 0
            if np.isnan(SN.dm15_array[i]):
                SN.dm15_array[i] = 0
            if SN.age[i] != 0:
                SN.age[i] = SN.phase
            if SN.dm15_array[i] != 0:
                SN.dm15_array[i] = SN.dm15
            if np.isnan(SN.dm15_array[i]):
                SN.dm15_array[i] = 0
	    if SN.red_array[i] != 0:
		SN.red_array[i] = SN.redshift	    
	    if np.isnan(SN.red_array[i]):
		SN.red_array[i] = 0
	    if SN.vel[i] != 0:
		SN.vel[i] = SN.velocity
	    if np.isnan(SN.vel[i]):
		SN.vel[i] = 0

        #print SN.filename
    #Here we clean up the data we pulled
    #Some supernovae are missing important data, so we just get rid of them
    #This can take a very long time if you have more than 500 spectra
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
    SN_Array = [SN for SN in SN_Array if hasattr(SN, 'ivar')]
    SN_Array = [SN for SN in SN_Array if SN.phase != None]
    SN_Array = [SN for SN in SN_Array if SN.redshift != None]
    print len(SN_Array), "spectra remain"
    return SN_Array


#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

#This is the model for how scales should be applied, used in the find_scales function
def scale_func(vars, in_data, out_data):
    
    scale  = vars[0]
    model  = scale * in_data
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
        overlap   = temp_ivar * ivar
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

            #print "Scale factor = ", result

            scales = np.append(scales, np.array([float(result)]), axis = 0)

    return scales

#Scales the data using the factors found before
#If a scale of zero is found, the spectrum's variance becomes zero so it just doesn't count.
badfiles = []
def scale_data(SN_Array, scales):
    print "Scaling..."
    for i in range(len(scales)):
        if scales[i] != 0:
            SN_Array[i].flux *= np.abs(scales[i])
            SN_Array[i].ivar /= (scales[i])**2
            #print "Scaled at factor ", scales[i]
        else:
            SN_Array[i].ivar = np.zeros(len(SN_Array[i].ivar))
    return SN_Array

#averages with weights of the inverse variances in the spectra
def average(SN_Array, template, medmean):
        print "Averaging..."
        #print fluxes, errors
        fluxes = []
        ivars  = []
        reds   = []
        phases = []
        ages   = []
        vels   = []
        dm15s  = []
        for SN in SN_Array:
            if len(fluxes) == 0:
                fluxes = np.array([SN.flux])
                ivars  = np.array([SN.ivar])
                reds   = np.array([SN.red_array])
                phases = np.array([SN.phase])
                ages   = np.array([SN.age])
                vels   = np.array([SN.vel])
                dm15s  = np.array([SN.dm15_array])
            else:
                try:
                    fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
                    ivars  = np.append(ivars, np.array([SN.ivar]), axis=0)
                    reds   = np.append(reds, np.array([SN.red_array]), axis = 0)
                    phases = np.append(phases, np.array([SN.phase]), axis = 0)
                    ages   = np.append(ages, np.array([SN.age]), axis = 0)
                    vels   = np.append(vels, np.array([SN.vel]), axis = 0)
                    dm15s  = np.append(dm15s, np.array([SN.dm15_array]), axis = 0)
                except ValueError:
                    print "This should never happen!"
        #Make flux/ivar mask so we can average for pixels where everything has 0 ivar
        flux_mask = np.zeros(len(fluxes[0,:]))
        ivar_mask = np.zeros(len(fluxes[0,:]))
        have_data = np.where(np.sum(ivars, axis = 0)>0)
        no_data   = np.where(np.sum(ivars, axis = 0)==0)
        ivar_mask[no_data] = 1

        #Add in flux/ivar mask
        fluxes = np.append(fluxes, np.array([flux_mask]), axis=0)
        ivars  = np.append(ivars, np.array([ivar_mask]), axis=0)
        reds   = np.append(reds, np.array([flux_mask]), axis=0)
        #phases = np.append(phases, np.array([flux_mask]), axis=0)
        ages   = np.append(ages, np.array([flux_mask]), axis=0)
        vels   = np.append(vels, np.array([flux_mask]), axis=0)
        dm15s  = np.append(dm15s, np.array([flux_mask]), axis=0)

        #The way this was done before was actually creating a single value array...
        #This keeps them intact correctly.
        reds      = [red for red in reds if red != None]
        phases    = [phase for phase in phases if phase != None]
        ages      = [age for age in ages if age != None]
        #vels      = [vel for vel in vels if vel != None]
        #vels      = [vel for vel in vels if vel != -99.0]
        dm15s     = [dm15 for dm15 in dm15s if dm15 != None]
        dm15_ivars = np.array(ivars)
############
# This is supposed to get rid of lists where there are no dm15 values
# The same logic will eventually be applied to phase and other single-value parameters that may be missing from some spectra
# It...doesn't work.
	"""
        for i in range(len(dm15s)):
            #print dm15s[i]
            if sum(dm15s[i]) == 0:
                #print "all zeros here"
                np.delete(dm15s, i, 1)
                np.delete(dm15_ivars, i, 1)
		
	#uh..this crashes my computer. don't use this bit.
	dm15s = dm15s[(dm15s != 0)-1]
	dm15_ivars = dm15_ivars[(dm15s != 0)-1]
	print dm15_ivars
	"""
	if len(dm15_ivars)==len(ivars):
            print "No rows were removed :("
############        
        if medmean == 1:
            template.flux  = np.average(fluxes, weights=ivars, axis=0)
            #template.phase = np.average(phases, weights=ivars, axis=0)
            template.age   = np.average(ages, weights=ivars, axis=0)
            template.vel   = np.average(vels, weights=ivars, axis=0)
            template.dm15  = np.average(dm15s, weights=dm15_ivars, axis=0)
	    template.red_array = np.average(reds, weights = ivars, axis=0)
        if medmean == 2:
            template.flux  = np.median(fluxes, axis=0)
            template.phase = np.median(phases, axis=0)
            template.age   = np.median(ages, axis=0)
            template.vel   = np.median(vels, axis=0)
            template.dm15  = np.median(dm15s, axis=0)
	    template.red_array = np.median(reds, axis=0)
        template.ivar = 1/np.sum(ivars, axis=0)
        #try:
        #    template.redshift = sum(reds)/len(reds)
        #except ZeroDivisionError:
        #    template.redshift = "No redshift data"
        try:
            template.phase = sum(phases)/len(phases)
        except ZeroDivisionError:
            template.phase = "No phase data"
        try:
            template.velocity = sum(vels)/len(vels)
        except ZeroDivisionError:
            template.velocity = "No velocity data"
        template.ivar[no_data] = 0
        template.name = "Composite Spectrum"
        return template



def main(Full_query, showplot = 0, medmean = 1, save_file = 'y'):
    SN_Array = []
    
    #Accept SQL query as input and then grab what we need
    print "SQL Query:", Full_query
    sql_input = Full_query

    SN_Array = grab(sql_input, Full_query)

    opt = str(raw_input("Do you want to do bootstrapping? (y/n) "))
    if opt == 'n':
	#finds the longest SN we have for our initial template
        lengths = []
        for SN in SN_Array:
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        try:
            composite = temp[0]
        except IndexError:
            print "No spectra found"
            exit()
        
    
        #scales data, makes a composite, and splices in non-overlapping data
        #Here is where we set our wavelength range for the final plot
        wmin    = 4000
        wmax    = 7500
        wavemin = composite.minwave
        wavemax = composite.maxwave
    
        #finds range of useable data
        good     = np.where(len(np.where(((wavemin <= wmin) & (wavemax >= wmax)) > 100)))
        template = supernova()
        template = SN_Array[good[0]]
        template = composite
        
        #Starts our main loop
        i = 0
        n_start = 0
        n_end   = 1
        scales  = []
        while (n_start != n_end):
            n_start = len([x for x in scales if x>0])
            scales   = []       
            scales   = find_scales(SN_Array, template.flux, template.ivar)
            n_scale  = len([x for x in scales if x>0])
            SN_Array = scale_data(SN_Array, scales)
            template = average(SN_Array, template, medmean)
            n_end    = n_scale
            n_start  = n_end
            
            
        print "Done."
        print "Average redshift =", template.redshift
        print "Average phase =", template.phase
        print "Average velocity =", template.velocity
        #This next line creates a unique filename for each run based on the sample set used
	#### file_name.py needs to be adjusted
        #f_name = "../plots/" + file_name.make_name(SN_Array)
	f_name = "../plots/" + "Test_composite"
        template.savedname = f_name + '.dat'
        lowindex  = np.where(template.wavelength == find_nearest(template.wavelength, wmin))
        highindex = np.where(template.wavelength == find_nearest(template.wavelength, wmax))
        
        #This plots the individual composite just so you can see how it 
        if int(showplot) == 1:
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.flux[lowindex[0]:highindex[0]])
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.ivar[lowindex[0]:highindex[0]])
        
            #This saves it, if you want to.
            plt.savefig('../plots/' + f_name + '.png')
            plt.show()
        #Either writes data to file, or returns it to user
        #This part is still in progress
        table = Table([template.wavelength, template.flux, template.ivar, template.age, template.vel, template.dm15], names = ('Wavelength', 'Flux', 'Variance', 'Age', 'Velocity', 'Dm_15'))
        if save_file=='y':
                    table.write(template.savedname,format='ascii')
                    return template
        else:
                    return template
    if opt == 'y':
        tries = int(raw_input("Enter number of bootstraps: "))  # Number of bootstraps.
        
        boot_flux = []
        boot_flux_unscaled = [0] * tries
        
        
        for j in range(tries):
            SN_Array2 = bootstrap.main(SN_Array)
            #finds the longest SN we have for our initial template
            lengths = []
            for SN in SN_Array2:
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
            temp = [SN for SN in SN_Array2 if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
            try:
                composite = temp[0]
            except IndexError:
                print "No spectra found"
                exit()
        
            #scales data, makes a composite, and splices in non-overlapping data
            #Here is where we set our wavelength range for the final plot
            wmin    = 4000
            wmax    = 7500
            wavemin = composite.minwave
            wavemax = composite.maxwave
        
            #finds range of useable data
            good     = np.where(len(np.where(((wavemin <= wmin) & (wavemax >= wmax)) > 100)))
            template = supernova()
            template = SN_Array2[good[0]]
            template = composite
        
            #Starts our main loop
            i = 0
            n_start = 0
            n_end   = 1
            scales  = []
            while (n_start != n_end):
                n_start = len([x for x in scales if x>0])
                scales   = []
                scales   = find_scales(SN_Array2, template.flux, template.ivar)
                n_scale  = len([x for x in scales if x>0])
                SN_Array2 = scale_data(SN_Array2, scales)
                template = average(SN_Array2, template, medmean)
                n_end    = n_scale
                n_start  = n_end
            
            
            boot_flux_unscaled[j] = template.flux[np.where(template.wavelength == wmin)[0]:np.where(template.wavelength == wmax)[0]]
        
            boot_flux.append(np.divide(boot_flux_unscaled[j], np.median(boot_flux_unscaled[j])))
                        
                
        
        
        ### 16th and 84th percentile of the spectrum (for scatter plot)
        percentile = erf(1/np.sqrt(2.))

        low_pc = 0.5 - percentile / 2.
        up_pc = 0.5 + percentile / 2.

        ### The 16th and 84th percentile index
        low_ind = np.round(tries * low_pc).astype(int)
        up_ind = np.round(tries * up_pc).astype(int)

        ### Sort the fluxes in each wavelength, and put the 16th and 84th percentile fluxes into two arrays
        median = np.median(boot_flux, axis = 0)     ### Median of the spectrum (for scaling)
        low_arr = np.divide(np.sort(boot_flux, axis = 0)[low_ind - 1], median)
        up_arr = np.divide(np.sort(boot_flux, axis = 0)[up_ind - 1], median)

            #for j in range(np.sort(boot_flux, axis = 0)[low_ind - 1].size):
            #print j, median[j], aaaa[j], low_arr[j], up_arr[j]
        
        lowindex  = np.where(template.wavelength == find_nearest(template.wavelength, wmin))
        highindex = np.where(template.wavelength == find_nearest(template.wavelength, wmax))
                        
            #print low_arr[lowindex[0]:highindex[0]]
                
        plt.plot(template.wavelength[lowindex[0]:highindex[0]], low_arr)
        plt.plot(template.wavelength[lowindex[0]:highindex[0]], up_arr)
        print low_arr
        print up_arr
        
        minflux = np.min(low_arr) * 0.9
        maxflux = np.max(up_arr) * 1.1
        plt.ylim((minflux, maxflux))
        plt.show()
        plt.close()
                
                #print median[lowindex[0]:highindex[0]]
        f_name = "../plots/" + file_name.make_name(SN_Array2)
        template.savedname = f_name + '.dat'
        #This plots the individual composite just so you can see how it
        if int(showplot) == 1:
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.flux[lowindex[0]:highindex[0]])
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.ivar[lowindex[0]:highindex[0]])
                
            #This saves it, if you want to.
            plt.savefig('../plots/' + f_name + '.png')
            plt.show()
        #Either writes data to file, or returns it to user
        #This part is still in progress
	print template.phase
        table = Table([template.wavelength, template.flux, template.ivar, template.age, template.dm15s], names = ('Wavelength', 'Flux', 'Variance', 'Age', 'Dm_15s'))
        if save_file=='y':
            table.write(template.savedname,format='ascii')
            return template
        else:
            return template
if __name__ == "__main__":
    main()
