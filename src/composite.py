"""
Spectra composite program
Authors: Sam, Yixian, Aaron
"""
#example: python Run.py 2 "SELECT * FROM Supernovae WHERE Carbon = 'A' AND Phase Between -6 and -4 AND Dm15 Between 0 and 2"
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
import bootstrap
import time
from lmfit import minimize, Parameters
import scipy.optimize as opt
import copy

np.set_printoptions(threshold=np.nan)
mn.patch()

#Sets up some lists for later
SN_Array = []
full_array = []
compare_spectrum = []

class supernova(object):
    """Attributes can be added"""

#Connect to database
#Make sure your file is in this location
con = sq3.connect('../data/SNe.db')
cur = con.cursor()

#Pulls in all columns from the database for the selected query
def grab(sql_input):
    print "Collecting data..."
    SN_Array = []
    cur.execute(sql_input)
    for row in cur:
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
        for i in range(len(SN_Array)-1):
            if SN_Array[i].name == SN_Array[i-1].name:
                if abs(SN_Array[i].phase) < abs(SN_Array[i-1].phase):
                    del SN_Array[i-1]
    print len(SN_Array), "spectra found"

    #Within the interpolated spectra there are a lot of 'NaN' values
    #Now they become zeros so things work right
    for SN in SN_Array:
        SN.phase_array = np.array(SN.flux)
        SN.dm15_array  = np.array(SN.flux)
        SN.red_array   = np.array(SN.flux)
        SN.vel         = np.array(SN.flux)
        for i in range(len(SN.flux)):
            #Check for NaN
            if np.isnan(SN.flux[i]):
                SN.flux[i]         = 0
                SN.ivar[i]         = 0
                SN.phase_array[i]  = 0
                SN.dm15_array[i]   = 0
                SN.red_array[i]    = 0
                SN.vel[i]          = 0
            
            #Set nonzero values to correct ones
            if SN.phase_array[i] != 0:
                if SN.phase != None:
                    SN.phase_array[i] = SN.phase
                else:
                    SN.phase_array[i] = 0
            if SN.dm15_array[i] != 0:
                if SN.dm15 != None:
                    SN.dm15_array[i] = SN.dm15
                else:
                    SN.dm15_array[i] = 0
            if SN.red_array[i] != 0:
                if SN.redshift != None:
                    SN.red_array[i] = SN.redshift
                else:
                    SN.red_array[i] = 0
            if SN.vel[i] != 0:
                if SN.velocity != None:
                    SN.vel[i] = SN.velocity
                else:
                    SN.vel[i] = 0
                    
    print "Arrays cleaned"
    return SN_Array
    
    
#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

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
	elif len(flux) == 0:
	    return scales

        else:
            #Otherwise, fit things
            vars = [1.0]
            #Find the appropriate values for scaling
            good      = np.where(overlap > 0)
	    good      = np.array(good[0])
	    
            flux2     = np.array(flux[good[0]:good[-1]])
            ivar2     = ivar[good[0]:good[-1]]
            tempflux2 = np.array(temp_flux[good[0]:good[-1]])
            tempivar2 = temp_ivar[good[0]:good[-1]]
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
def scale_data(SN_Array, scales): #obselete
    print "Scaling..."
    print scales
    for i in range(len(scales)):
        if scales[i] != 0:
            SN_Array[i].flux *= np.abs(scales[i])
            SN_Array[i].ivar /= (scales[i])**2
            #print "Scaled at factor ", scales[i]
        else:
            SN_Array[i].ivar = np.zeros(len(SN_Array[i].ivar))
	#print SN_Array[i].flux
	if len(SN_Array[i].flux) == 0:
	    np.delete(SN_Array,i)
	else:
	    print len(SN_Array[i].wavelength), len(SN_Array[i].flux)
	    plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
    print "Displaying scaled spectra..."
    plt.show()
    return SN_Array


#########
#Here's the new scale function that I'm trying
#Based on some code that Brian sent me
#I think it works the way it's designed to at the moment
#But I think it might not be what we want.
#It fits the spectra to each other...
#but instead of scaling them while remaining the same shape,
#they all get fitted to a single curve
#So they end up identical...
#bad, right?
##########

def new_scale_func(SN_Array, template): #obselete
    #possibly use scipy.optimize.minimize 
    #write function that sums all residuals so minimize solves for coeffs
    print "Using fitting function..."
    for SN in SN_Array:
        data = np.array([template.flux, SN.flux])
        source  = data.T
        print len(source)
	if len(source) < 3000:
	    scales = []
	    return SN_Array, scales
        guess = 4
        datas = []
        for m in range(len(source)):
            comp2, cov2 = leastsq(residual, guess, args=(source[m]), full_output=False)
            datas.append(comp2[0])
        #print datas
        scales = [x for x in datas if x != 0]
        SN.flux = datas
        plt.plot(SN.wavelength, SN.flux)
    plt.show()
    return SN_Array, scales       
    
    
def optimize_scales(SN_Array, template):
    print "Minimizing Square Residuals..."
    for SN in SN_Array:
        SN.flux = SN.flux*1e14 #values too small otherwise?
    template.flux = template.flux*1e14
    guess = np.ones(len(SN_Array)) # this must be ones so the template is multiplied by 1
    scales = opt.minimize(sq_residuals, guess, args = (SN_Array, template)).x
    print scales
    for i in range(len(SN_Array)):
        SN_Array[i].flux = 1e-14*scales[i]*SN_Array[i].flux
        SN_Array[i].ivar /= (scales[i])**2
        plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
    template.flux = 1e-14*template.flux
    plt.plot(template.wavelength, template.flux,'k', linewidth = 4)
    plt.show()
    return SN_Array, scales
        
        
def sq_residuals(s, SN_Array, comp):
    all_sums = []
    for i in range(len(SN_Array)):
#        if SN_Array[i] != comp:
        SN = SN_Array[i]
        non_zero_data = np.where(SN.flux != 0)
        non_zero_data = np.array(non_zero_data[0])
        
        non_zero_data_comp = np.where(comp.flux != 0)
        non_zero_data_comp = np.array(non_zero_data_comp[0])
        
        if len(set(non_zero_data).intersection(non_zero_data_comp)) != 0:
            temp_sum = 0
            pos1 = non_zero_data[0]
            pos2 = non_zero_data[-1]
            temp_flux = s[i]*SN.flux
            res = comp.flux[pos1:pos2] - temp_flux[pos1:pos2]
            sq_res = res**2
            temp_sum = np.sum(sq_res)
            all_sums.append(temp_sum)
    
    tot = np.sum(all_sums)
    return tot
    
def residual(comp, data):
    return comp-data

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
                ages   = np.array([SN.phase_array])
                vels   = np.array([SN.vel])
                dm15s  = np.array([SN.dm15_array])
            else:
                try:
                    fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
                    ivars  = np.append(ivars, np.array([SN.ivar]), axis=0)
                    reds   = np.append(reds, np.array([SN.red_array]), axis = 0)
                    phases = np.append(phases, np.array([SN.phase]), axis = 0)
                    ages   = np.append(ages, np.array([SN.phase_array]), axis = 0)
                    vels   = np.append(vels, np.array([SN.vel]), axis = 0)
                    dm15s  = np.append(dm15s, np.array([SN.dm15_array]), axis = 0)
                except ValueError:
                    print "This should never happen!"

        
        #Adding masks for every parameter for consistency and zero compensation
        flux_mask = np.zeros(len(fluxes[0,:]))
        ivar_mask = np.zeros(len(fluxes[0,:]))
        dm15_mask = np.zeros(len(dm15s[0,:]))
        red_mask  = np.zeros(len(reds[0,:]))
        
        have_data = np.where(np.sum(ivars, axis = 0)>0)
        no_data   = np.where(np.sum(ivars, axis = 0)==0)
        no_dm15   = np.where(np.sum(dm15s, axis = 0)==0)
        no_reds   = np.where(np.sum(reds, axis = 0)==0)
        
        ivar_mask[no_data] = 1
        dm15_mask[no_dm15] = 1
        
        #Right now all of our spectra have redshift data, so a mask is unnecessary
        #One day that might change?
        red_mask[:]  = 1
        
        dm15_ivars = np.array(ivars)
        red_ivars  = np.array(ivars)
        
        #Add in flux/ivar mask
        fluxes = np.append(fluxes, np.array([flux_mask]), axis=0)
        ivars  = np.append(ivars, np.array([ivar_mask]), axis=0)
        reds   = np.append(reds, np.array([flux_mask]), axis=0)
        ages   = np.append(ages, np.array([flux_mask]), axis=0)
        vels   = np.append(vels, np.array([flux_mask]), axis=0)
        dm15s  = np.append(dm15s, np.array([dm15_mask]), axis=0)
        dm15_ivars = np.append(dm15_ivars, np.array([dm15_mask]), axis=0)
        red_ivars  = np.append(red_ivars, np.array([red_mask]), axis=0)

        for i in range(len(dm15s)):
            if np.all(dm15s[i]) == 0:
                np.delete(dm15s, i)
                np.delete(dm15_ivars, i)

        if medmean == 1: 
            template.flux  = np.average(fluxes, weights=ivars, axis=0)
            template.phase_array   = np.average(ages, weights=ivars, axis=0)
            template.vel   = np.average(vels, weights=ivars, axis=0)
            template.dm15  = np.average(dm15s, weights=dm15_ivars, axis=0)
            template.red_array = np.average(np.array(reds), weights = red_ivars, axis=0)
        if medmean == 2:
            template.flux  = np.median(fluxes, axis=0)
            template.phase_array   = np.median(ages, axis=0)
            template.vel   = np.median(vels, axis=0)
            template.dm15  = np.median(dm15s, axis=0)
            template.red_array = np.median(reds, axis=0)
            
        template.ivar = 1/np.sum(ivars, axis=0)
        template.ivar[no_data] = 0
        template.name = "Composite Spectrum"
        return template

def do_things(SN_Array, template, scales, medmean):
    n_start = len([x for x in scales if x>0])
    scales   = []
    n_end = n_start
#    scales   = find_scales(SN_Array, template.flux, template.ivar)
#    SN_Array = scale_data(SN_Array, scales)
    SN_Array, scales = optimize_scales(SN_Array, template)
#   SN_Array, scales = new_scale_func(SN_Array, template)
    template = average(SN_Array, template, medmean)
    return n_start, n_end, SN_Array, template, scales
        
def bootstrapping (SN_Array, trials, scales, og_template):
    strap_matrix = np.random.random_sample((trials, len(SN_Array)))
    strap_matrix *= len(SN_Array)-1
    strap_matrix = np.round(strap_matrix)
    strap_matrix = strap_matrix.astype(int)
    
    boot_arr = []
    boots = []
    for i in range(len(strap_matrix)):
        boot_arr.append([])
        for j in range(len(strap_matrix[i])):
            boot_arr[i].append(copy.copy(SN_Array[strap_matrix[i,j]]))

    for p in range(len(boot_arr)):
        lengths = []
        for SN in boot_arr[p]:
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        boot_temp = [SN for SN in boot_arr[p] if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        boot_temp = copy.copy(boot_temp[0])
    
        for x in range(3):
            n_start, n_end, SN_Array, template, scales = do_things(boot_arr[p], boot_temp, scales, 1)
        boots.append(copy.copy(template))
    
    print "scaling boots..."
    optimize_scales(boots, og_template) 
    boot_wave = []
    boot_flux = []
    for SN in boots:
        boot_wave.append(SN.wavelength)
        boot_flux.append(SN.flux)
        
    plt.plot(boot_wave, boot_flux)
    plt.show()
    

def main(Full_query, showplot = 0, medmean = 1, opt = 'n', save_file = 'n'):
    SN_Array = []

    #Accept SQL query as input and then grab what we need
    print "SQL Query:", Full_query
    SN_Array = grab(Full_query)
    
    if opt == 'y':
        bootstrap.main(SN_Array)

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

        #Here is where we set our wavelength range for the final plot
        wmin    = 4000
        wmax    = 7500
        wavemin = composite.minwave
        wavemax = composite.maxwave

        #finds range of useable data
        good     = np.where(len(np.where(((wavemin <= wmin) & (wavemax >= wmax)) > 100)))
        template = supernova()
        template = copy.copy(composite)

        #Starts our main loop
        #scales data, makes a composite, and splices in non-overlapping data
        i = 0
        n_start = 0
        n_end   = 1
        scales  = []
        
        print "Creating composite..."
#        while (n_start != n_end):
        for i in range(3):
            n_start, n_end, SN_Array, template, scales = do_things(SN_Array, template, scales, medmean)
        n_start, n_end, SN_Array, template, scales = do_things(SN_Array, template, scales, medmean)
        print "Done."

        scales  = []
#        print "Bootstrapping"
#        bootstrapping(SN_Array, 10, scales, template)
        #This next line creates a unique filename for each run based on the sample set used
        #### file_name.py needs to be adjusted
        #f_name = "../plots/" + file_name.make_name(SN_Array)
        f_name = "../plots/" + "Test_composite" + (time.strftime("%H,%M,%S"))
        template.savedname = f_name + '.dat'
        lowindex  = np.where(template.wavelength == find_nearest(template.wavelength, wmin))
        highindex = np.where(template.wavelength == find_nearest(template.wavelength, wmax))

        #This plots the individual composite just so you can see how it looks
        #Also it gets saved. Comment that line out if you don't want it to.
        if int(showplot) == 1:
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.flux[lowindex[0]:highindex[0]])
            plt.plot(template.wavelength[lowindex[0]:highindex[0]], template.ivar[lowindex[0]:highindex[0]])
#            plt.savefig('../plots/' + f_name + '.png')
            plt.show()
        table = Table([template.wavelength, template.flux, template.ivar, template.phase_array, template.vel, template.dm15, template.red_array], names = ('Wavelength', 'Flux', 'Variance', 'Age', 'Velocity', 'Dm_15', 'Redshift'))
        if save_file == 'y':
            table.write(template.savedname, format='ascii')
        return table

if __name__ == "__main__":
    main()
