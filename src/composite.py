"""
Spectra composite program
Authors: Sam, Yixian, Aaron
"""
#example: python Run.py 2 "SELECT * FROM Supernovae WHERE Carbon = 'A' AND Phase Between -6 and -4 AND Dm15 Between 0 and 2"
#msgpack_python version 0.4.6
#msgpack_numpy version 0.3.5
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
import time
from lmfit import minimize, Parameters
import scipy.optimize as opt
import copy
from collections import Counter

from specutils import extinction as ex
import test_dered
import prep
from astropy import units as u
from specutils import Spectrum1D
import magnitudes as mg

np.set_printoptions(threshold=np.nan)
mn.patch()

#Sets up some lists for later
SN_Array = []
full_array = []
compare_spectrum = []


class supernova(object):
    """Contains all spectral data provided by the associated file.
       Attributes can be added
    """

#Connect to database
#Make sure your file is in this location

#con = sq3.connect('../data/SNe.db')
# con = sq3.connect('../data/SNe_2.db')
con = sq3.connect('../data/SNe_3.db')
cur = con.cursor()


def grab(sql_input):
    """Pulls in all columns from the database for the selected query. 
       Replaces all NaN values with 0. Returns the array of supernova objects 
       with the newly added attributes.
    """
    print "Collecting data..."
    SN_Array = []
    multi_epoch = raw_input("Include multiple epochs? (y/n): ")
    if multi_epoch == 'y':
        multi_epoch = True
    else:
        multi_epoch = False

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
        SN.resid     = row[15]
        interp       = msg.unpackb(row[16])
        SN.interp    = interp
        phot         = msg.unpackb(row[17])
        SN.phot      = phot
        try:
            SN.wavelength = SN.interp[0,:]
            SN.flux       = SN.interp[1,:]
            SN.ivar       = SN.interp[2,:]
        except TypeError:
            continue
        full_array.append(SN)
        SN_Array.append(SN)

        # for i in range(len(SN_Array-1)): 
        #     if SN_Array[i].name == SN_Array[i-1].name and not multi_epoch:
        #         if abs(SN_Array[i].phase) < abs(SN_Array[i-1].phase): # closest to maximum
        #         # if abs(SN_Array[i].SNR) > abs(SN_Array[i-1].SNR): # best signal to noise
        #             del SN_Array[i-1]
        if not multi_epoch:
            unique_events = []
            new_SN_Array = []
            for i in range(len(SN_Array)): 
                if SN_Array[i].name not in unique_events:
                    unique_events.append(SN_Array[i].name)
            for i in range(len(unique_events)):
                events = []
                for SN in SN_Array:
                    if SN.name == unique_events[i]:
                        events.append(SN)
                min_phase = events[0]
                for e in events:
                    if abs(e.phase) < abs(min_phase.phase):
                        min_phase = e
                new_SN_Array.append(min_phase)
            SN_Array = new_SN_Array

    print len(SN_Array), "spectra found"

    # print "Creating event file..."
    # mg.generate_event_list(SN_Array)
    # print "Event file done."
    
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
#        print SN.ivar
        non_zero_data = np.where(SN.flux != 0)
        non_zero_data = np.array(non_zero_data[0])
        if len(non_zero_data) > 0:
            SN.x1 = non_zero_data[0]
            SN.x2 = non_zero_data[-1]
        else:
            SN.x1 = 0.
            SN.x2 = 0.
                    
    print "Arrays cleaned"
    return SN_Array

    
def spectra_per_bin(SN_Array):
    """Counts the number of spectra contributing to the composite at any given 
       wavelength. Returns array to be plotted over the wavelength range.
    """
    spec_per_bin = []
    
    for i in range(len(SN_Array[0].flux)):
        count = 0
        for SN in SN_Array:
            if SN.flux[i] != 0 and SN.ivar[i] != 0:   
                count += 1
        spec_per_bin.append(count)
    
    return spec_per_bin
            
            
def optimize_scales(SN_Array, template, initial):
    """Scales each unique supernova in SN_Array by minimizing the square residuals
       between the supernova flux and the template flux. This also works for bootstrap
       arrays (can contain repeated spectra) because the objects in SN_Array are not 
       copies. Returns scaled SN_Array and the scales that were used.
    """
    scales = []
    unique_arr = list(set(SN_Array))
    guess = 1.0
    for uSN in unique_arr:
#        guess = template.flux[2000]/uSN.flux[2000]
        if uSN.filename != template.filename:
            u = opt.minimize(sq_residuals, guess, args = (uSN, template, initial), 
                             method = 'Nelder-Mead').x
            scales.append(u)
        else:
            scales.append(1.0)
        
    for i in range(len(unique_arr)):
        unique_arr[i].flux = scales[i]*unique_arr[i].flux
        unique_arr[i].ivar /= (scales[i])**2
        
    return SN_Array, scales
    
    
def sq_residuals(s,SN,comp, initial):
    """Calculates the sum of the square residuals between two supernova flux 
       arrays usinig a given scale s. Returns the sum.
    """
#    print s
    if SN.x1 <= comp.x1 and SN.x2 >= comp.x2:
        pos1 = comp.x1
        pos2 = comp.x2
    elif SN.x1 >= comp.x1 and SN.x2 <= comp.x2:
        pos1 = SN.x1
        pos2 = SN.x2
    elif SN.x1 >= comp.x1 and SN.x1 <= comp.x2 and SN.x2 >= comp.x2:
        pos1 = SN.x1
        pos2 = comp.x2
    elif SN.x1 <= comp.x1 and SN.x2 >= comp.x1 and SN.x2 <= comp.x2:
        pos1 = comp.x1
        pos2 = SN.x2
    else: 
        #residual of first index will always be zero (won't be scaled)
        print "no overlap "
        pos1 = 0
        pos2 = 0
    temp_flux = s*SN.flux
    res = temp_flux[pos1:pos2] - comp.flux[pos1:pos2]
    sq_res = res*res
    if initial:
        return np.sum(sq_res)
    else:
        temp_ivar = SN.ivar/(s*s)
        w_res = temp_ivar[pos1:pos2]*sq_res
        return np.sum(w_res)
    
def mask(SN_Array, boot):
    """Creates data structures to contain relevant data for the task needed 
       (creating the composite or bootstrapping). Applies masks the these data 
       for consistency and zero compensation. Returns the masks and data structures.
    """
    #create 2D arrays of all available data
    fluxes = []
    ivars  = []
    reds   = []
    phases = []
    ages   = []
    vels   = []
    dm15s  = []
    dm15_ivars = []
    red_ivars  = []
    dm15_mask  = []
    red_mask   = []
    for SN in SN_Array:
        if len(fluxes) == 0:
            fluxes = np.array([SN.flux])
            ivars  = np.array([SN.ivar])
            if not boot:
                reds   = np.array([SN.red_array])
                phases = np.array([SN.phase])
                ages   = np.array([SN.phase_array])
                vels   = np.array([SN.vel])
                dm15s  = np.array([SN.dm15_array])
        else:
            try:
                fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
                ivars  = np.append(ivars, np.array([SN.ivar]), axis=0)
                if not boot:
                    reds   = np.append(reds, np.array([SN.red_array]), axis = 0)
                    phases = np.append(phases, np.array([SN.phase]), axis = 0)
                    ages   = np.append(ages, np.array([SN.phase_array]), axis = 0)
                    vels   = np.append(vels, np.array([SN.vel]), axis = 0)
                    dm15s  = np.append(dm15s, np.array([SN.dm15_array]), axis = 0)
            except ValueError:
                print "This should never happen!"

    #Adding masks for every parameter 
    flux_mask = np.zeros(len(fluxes[0,:]))
    ivar_mask = np.zeros(len(fluxes[0,:]))
    if not boot:
        dm15_mask = np.zeros(len(dm15s[0,:]))
        red_mask  = np.zeros(len(reds[0,:]))
    
    have_data = np.where(np.sum(ivars, axis = 0)>0)
    no_data   = np.where(np.sum(ivars, axis = 0)==0)
    if not boot:
        no_dm15   = np.where(np.sum(dm15s, axis = 0)==0)
        no_reds   = np.where(np.sum(reds, axis = 0)==0)
        dm15_mask[no_dm15] = 1
    
    ivar_mask[no_data] = 10e6
    
    #Right now all of our spectra have redshift data, so a mask is unnecessary
    #One day that might change
    if not boot:
        red_mask[:]  = 1
        
        dm15_ivars = np.array(ivars)
        red_ivars  = np.array(ivars)
    
    #Add in flux/ivar mask
    fluxes = np.append(fluxes, np.array([flux_mask]), axis=0)
    ivars  = np.append(ivars, np.array([ivar_mask]), axis=0)
    if not boot:
        reds   = np.append(reds, np.array([flux_mask]), axis=0)
        ages   = np.append(ages, np.array([flux_mask]), axis=0)
        vels   = np.append(vels, np.array([flux_mask]), axis=0)
        dm15s  = np.append(dm15s, np.array([dm15_mask]), axis=0)
        dm15_ivars = np.append(dm15_ivars, np.array([dm15_mask]), axis=0)
        red_ivars  = np.append(red_ivars, np.array([red_mask]), axis=0)
        
        for i in range(len(dm15_ivars)):
            for j in range(len(dm15_ivars[i])):
                if dm15_ivars[i,j] == 0.0:
                    dm15_ivars[i,j] = 1.0
    
    return (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, 
            dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
                

def average(SN_Array, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask):
    """Modifies the template supernova to be the inverse variance weighted average
       of the scaled data. Returns the new template supernova. 
    """
    for i in range(len(SN_Array)):
        fluxes[i] = SN_Array[i].flux
        ivars[i] = SN_Array[i].ivar            
    if medmean == 1: 
        template.flux  = np.average(fluxes, weights=ivars, axis=0)
        if not boot:
            template.phase_array   = np.average(ages, weights=ivars, axis=0)
            template.vel   = np.average(vels, weights=ivars, axis=0)
            template.dm15  = np.average(dm15s, weights=dm15_ivars, axis=0)
            template.red_array = np.average(np.array(reds), weights = red_ivars, axis=0)
    if medmean == 2:
        template.flux  = np.median(fluxes, axis=0)
        if not boot:
            template.phase_array   = np.median(ages, axis=0)
            template.vel   = np.median(vels, axis=0)
            template.dm15  = np.median(dm15s, axis=0)
            template.red_array = np.median(reds, axis=0)
    
    #finds and stores the variance data of the template
    no_data   = np.where(np.sum(ivars, axis = 0)==0)
    template.ivar = 1/np.sum(ivars, axis=0) #change to variance
    template.ivar[no_data] = 0
    template.name = "Composite Spectrum"
    return template

        
def bootstrapping (SN_Array, samples, scales, og_template, iters):
    """Creates a matrix of random sets of supernovae from the original sample 
       with the same size as the original sample. The number of samples is 
       defined by the user. Then creates and plots the composite spectrum for 
       each of these sets. These data are used to contruct a confidence 
       interval for the original sample. Returns flux arrays corresponding top 
       the upper and lower intervals.
    """
    strap_matrix = np.random.random_sample((samples, len(SN_Array)))
    strap_matrix *= len(SN_Array)
    strap_matrix = strap_matrix.astype(int)   
    boot_arr = []
    boots = []
    boot = True
    
    cpy_array = []
    for SN in SN_Array:
        cpy_array.append(copy.copy(SN))
        
    
    for i in range(len(strap_matrix)):
        boot_arr.append([])
        for j in range(len(strap_matrix[i])):
            boot_arr[i].append(SN_Array[strap_matrix[i,j]])
            

    for p in range(len(boot_arr)):
        print p
        lengths = []
        for SN in boot_arr[p]:
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        boot_temp = [SN for SN in boot_arr[p] if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        boot_temp = copy.copy(boot_temp[0])
        
        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
         flux_mask, ivar_mask, dm15_mask, red_mask) = mask(boot_arr[p], boot)
        for x in range(iters):
            SN_Array, scales = optimize_scales(boot_arr[p], boot_temp, False)
            template = average(boot_arr[p], boot_temp, 1, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                               reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
        boots.append(copy.copy(template))

    
    print "scaling boots..."
    optimize_scales(boots, og_template, False)
    print "plotting..."
    for SN in boots:
        plt.plot(SN.wavelength, SN.flux, 'g')
    plt.plot(og_template.wavelength,og_template.flux, 'k', linewidth = 4)
    plt.show()
    
    print "computing confidence intervals..."
    resid = []
    percentile = erf(1/np.sqrt(2.))
    low_pc = 0.5 - percentile*0.5
    up_pc = 0.5 + percentile*0.5

    for SN in boots:
#        temp = SN.flux - og_template.flux
        resid.append(SN.flux - og_template.flux)
#        non_olap = np.where(og_template.flux + temp == 0.0)
#        temp = np.delete(temp, non_olap)
#        resid.append(temp)
        
    resid_trans = np.transpose(resid)
    resid_sort = np.sort(resid_trans)
    arr = []
    new_resid = []
    for i in range(len(resid_sort)):
        non_olap = []
        for j in range (len(resid_sort[i])):
            if resid_sort[i][j] + og_template.flux[i] == 0.0 and og_template.flux[i] != 0.0:
                non_olap.append(j)
        new_resid.append(np.delete(resid_sort[i],non_olap))
                
    for elem in resid_sort:
        for i in range (len(elem)):
            arr.append(elem[i])

    plt.hist(arr, 100)
    plt.show()
    
    low_arr = []
    up_arr = []
    for i in range(len(new_resid)):
        low_ind = np.round((len(new_resid[i])-1) * low_pc).astype(int)
        up_ind = np.round((len(new_resid[i])-1) * up_pc).astype(int)
        low_arr.append(og_template.flux[i] + new_resid[i][low_ind])
        up_arr.append(og_template.flux[i] + new_resid[i][up_ind])
        
#    for i in range(len(resid_sort)):
#        print len(resid_sort[i])
#        low_ind = np.round((len(resid_sort[i])-1) * low_pc).astype(int)
#        up_ind = np.round((len(resid_sort[i])-1) * up_pc).astype(int)
#        low_arr.append(og_template.flux[i] + resid_sort[i][low_ind])
#        up_arr.append(og_template.flux[i] + resid_sort[i][up_ind])
    
    plt.plot(og_template.wavelength, og_template.flux, 'k', linewidth = 4)
    plt.fill_between(og_template.wavelength, low_arr, up_arr, color = 'green')
    plt.show()
    
    return low_arr, up_arr        
    
def is_bad_data(SN, bad_files, reddened_spectra):
    for el in bad_files:
        if SN.filename == el:
            print 'bad file removed'
            return True
    for el in reddened_spectra:
        if SN.filename == el:
            print 'reddened spectrum removed'
            return True
    return False
    
    
def main(Full_query, showplot = 0, medmean = 1, opt = 'n', save_file = 'n'):
    """Main function. Finds supernovae that agree with user query, prompts user 
       on whether to bootstrap or just create a composite, then does so and stores 
       returns the relevant data for plotting in a table.
    """
    SN_Array = []

    #Accept SQL query as input and then grab what we need
    print "SQL Query:", Full_query
    SN_Array = grab(Full_query)

    #finds the longest SN we have for our initial template
    lengths = []
    bad_files =        ['sn2004ef-20040915.30-fast.flm', 'sn1994T-19940611.21-fast.flm',
                        'sn2006cj-20060523.33-fast.flm', 'sn1996ab-19960522.37-fast.flm', 
                        'sn1997bp-19970411.30-fast.flm', 'sn1995bd-19951225.27-fast.flm',
                        'sn1991t-19910418.flm'] 
                        #sn1991t-19910418.flm no si feature
                        #sn1996ab-19960522.37-fast.flm  has large negative values
                        #sn2004ef-20040915.30-fast.flm variance in database are all nan
                        #sn2006cj-20060523.33-fast.flm variance in database are all nan
                        #sn1995bd-19951225.27-fast.flm variance in database are all nan
                        #sn1997bp-19970411.30-fast.flm whole spectrum clearly shifted bluer
                        #sn2007af-20070314.44-fast.flm has very low variance which strongly biases results
                        #sn2001da-20010715.47-mmt.flm has a large wavelength range
    reddened_spectra = ['sn2003cg-20030331.21-fast.flm', 'sn2002cd-2n0020419.48-fast.flm',
                        'sn1996ai-19960621.23-fast.flm', 'sn1997dt-19971204.11-fast.flm',
                        'sn2006br-20060427.33-fast.flm', 'sn2003W-20030209.35-fast.flm',
                        'sn2004gs-20041216.49-fast.flm', 'sn2002fb-20020912.44-fast.flm',
                        'sn2005bo-20050418.21-fast.flm', 'sn1998de-19980801.41-fast.flm',
                        'sn1998bp-19980503.45-fast.flm', 'sn2005ke-20051125.30-fast.flm',
                        'sn2006cm-20060529.46-fast.flm', 'sn1997bp-19970409.29-fast.flm',
                        'sn2000cp-20000624.34-fast.flm', 'sn2005A-20050108.13-fast.flm',
                        'sn1995E-19950226.27-fast.flm',  'sn1999cl-19990612.17-fast.flm',
                        'sn2003cg-20030330.23-fast.flm', 'sn2006gj-20060920.43-fast.flm',
                        'sn2007ax-20070327.26-fast.flm', 'sn2006bz-20060506.16-fast.flm',
                        'sn2005A-20050110.11-ldss2.flm', 'sn1995E-19950228.23-fast.flm',
                        'sn1999cl-19990614.18-fast.flm', 'sn2006X-20060221.40-fast.flm',
                        'sn1999gd-19991208.52-fast.flm', 'sn2003cg-20030401.22-fast.flm',
                        'sn1996ai-19960620.15-mmt.flm',

                        #phase -15.0 to -12.0
                        # 'sn1995bd-19951223.34-fast.flm', 'sn2002bo-20020310.26-fast.flm',
                        # 'sn2002bo-20020311.23-fast.flm', 'SN07S_070131_b01_NTT_EM.dat',

                        #phase -11.99 to -9.0
                        'sn2002bo-20020311-ui-corrected.flm', 'sn2006gz-20060930.13-fast.flm',
                        'sn2007bm-20070423.23-fast.flm', 'sn2003W-20030129.35-mmt.flm',
                        'sn2006qo-20061201.436-ui.flm', 'sn1997dt-19971123.19-fast.flm',
                        'sn2005cf-20050601.385-ui.flm', 'sn2002dj-20020615.17-fast.flm',
                        'sn2007le-20071015.325-br-corrected.flm', 'sn2006cc-20060508.34-fast.flm',
                        'sn1997dt-19971124.09-fast.flm', 'sn2004bw-20040527.362-lowopt.flm',
                        'sn1999dq-19990905.45-fast.flm']
    
    good_SN_Array = [SN for SN in SN_Array if not is_bad_data(SN, bad_files, reddened_spectra)]
    SN_Array = good_SN_Array

    # good_SNs = []
    # for SN in SN_Array:
    #     plt.plot(SN.wavelength, SN.flux)
    #     plt.show()
    #     accept = raw_input("Accept? (y/n): ")
    #     if accept is 'y':
    #         print "Added"
    #         good_SNs.append(SN)

    # SN_Array = good_SNs
    # print len(SN_Array)
    
    # mags = np.sort(mg.ab_mags(SN_Array), axis = 0)
    # for i in range(len(mags)):
    #     print mags[i][0], mags[i][1] 


    for SN in SN_Array:
        # print SN.resid
        print SN.name, SN.filename
        host_reddened = prep.ReadExtin('../data/info_files/ryan_av.txt')
        old_wave = SN.wavelength*u.Angstrom        # wavelengths
        old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
        spec1d = Spectrum1D.from_array(old_wave, old_flux)
        # new_flux = test_dered.dered(sne, SN.name, spec1d.wavelength, spec1d.flux)
        new_flux = test_dered.host_correction(host_reddened, SN.name, old_wave, old_flux)
        SN.flux = new_flux.value
        lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

    temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
    try:
        composite = temp[0]
    except IndexError:
        print "No spectra found"
        exit()
                
    for i in range (len(SN_Array)):
        for j in range(len(SN.flux)):
            
#            SN_Array[i].ivar[j] = np.median(SN_Array[i].ivar[SN_Array[i].x1:SN_Array[i].x2])
            if j > SN_Array[i].x1 and j < SN_Array[i].x2:
                if j < SN_Array[i].x1 + 5:
                    SN_Array[i].ivar[j] = np.median(SN_Array[i].ivar[j:j+10])
                if j > SN_Array[i].x2 + 5:
                    SN_Array[i].ivar[j] = np.median(SN_Array[i].ivar[j-10:j])
                else:
                    SN_Array[i].ivar[j] = np.median(SN_Array[i].ivar[j-5:j+5])
            else:
                SN_Array[i].ivar[j] = 0.0
                
    spec_bin = spectra_per_bin(SN_Array)

    #finds range of useable data
    template = supernova()
    template = copy.copy(composite)

    #creates main composite
    i = 0
    scales  = []
    iters = 3
    iters_comp = 3
    boot = False
    
    cpy_array = []
    for SN in SN_Array:
        cpy_array.append(copy.copy(SN))
    
    # plt.figure(num = 2, dpi = 100, figsize = [30, 20], facecolor = 'w')
    for i in range(len(SN_Array)):
        plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
    plt.plot(template.wavelength, template.flux, 'k', linewidth = 4)
    plt.show()
        
    bootstrap = raw_input("Bootstrap? (y/n): ")
    print "Creating composite..."
    optimize_scales(SN_Array, template, True)
    (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
     flux_mask, ivar_mask, dm15_mask, red_mask) = mask(SN_Array, boot)
    
    for i in range(iters_comp):
        SN_Array, scales = optimize_scales(SN_Array, template, False)
        template = average(SN_Array, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                           reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
    print "Done."
    
    #plot composite with the scaled spectra
    # plt.figure(num = 2, dpi = 100, figsize = [30, 20], facecolor = 'w')
    if bootstrap is 'n':
        for i in range(len(SN_Array)):
            plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
        plt.plot(template.wavelength, template.flux, 'k', linewidth = 4)
        plt.show()
        low_conf = template.flux
        up_conf = template.flux
    #create bootstrap composites
    else:
        scales  = []
        print "Bootstrapping"
        samples = int (raw_input("# of samples:"))
        low_conf, up_conf = bootstrapping(SN_Array, samples, scales, template, iters)
    
    table = Table([template.wavelength, template.flux, template.ivar, template.phase_array, 
                   template.vel, template.dm15, template.red_array, low_conf, up_conf, spec_bin], 
                   names = ('Wavelength', 'Flux', 'Variance', 'Age', 'Velocity', 
                            'Dm_15', 'Redshift', 'Lower Confidence', 'Upper Confidence', 'Spectra Per Bin'))
    if save_file == 'y':
        table.write(template.savedname, format='ascii')
    return table

if __name__ == "__main__":
    main()
