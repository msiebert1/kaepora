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
from collections import Counter

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
con = sq3.connect('../data/SNe.db')
cur = con.cursor()


def grab(sql_input):
    """Pulls in all columns from the database for the selected query. 
       Replaces all NaN values with 0. Returns the array of supernova objects 
       with the newly added attributes.
    """
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
        
        non_zero_data = np.where(SN.flux != 0)
        non_zero_data = np.array(non_zero_data[0])
        SN.x1 = non_zero_data[0]
        SN.x2 = non_zero_data[-1]
                    
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
            if SN.flux[i] != 0:
                count += 1
        spec_per_bin.append(count)
    
    return spec_per_bin
            
            
def optimize_scales(SN_Array, template):
    """Scales each unique supernova in SN_Array by minimizing the square residuals
       between the supernova flux and the template flux. This also works for bootstrap
       arrays (can contain repeated spectra) because the objects in SN_Array are not 
       copies. Returns scaled SN_Array and the scales that were used.
    """
    scales = []
    guess = 1.0 
    unique_arr = list(set(SN_Array))
    for uSN in unique_arr:
        u = opt.minimize(sq_residuals, guess, args = (uSN, template)).x
        scales.append(u)
        
    for i in range(len(unique_arr)):
        unique_arr[i].flux = scales[i]*unique_arr[i].flux
        unique_arr[i].ivar /= (scales[i])**2
        
    return SN_Array, scales
    
def optimize_scales_old(SN_Array, template):
    """Scales each unique supernova in SN_Array by minimizing the square residuals
       between the supernova flux and the template flux. Only works for bootstrap 
       arrays if the array contains copies of each supernova object. Returns 
       scaled SN_Array and the scales that were used.
    """
    scales = []
    guess = 1.0 
    
    for SN in SN_Array:
        s = opt.minimize(sq_residuals, guess, args = (SN, template)).x
        scales.append(s)
        
    for i in range(len(SN_Array)):
        SN_Array[i].flux = scales[i]*SN_Array[i].flux
        SN_Array[i].ivar /= (scales[i])**2
        
    return SN_Array, scales
    
    
def sq_residuals(s,SN,comp):
    """Calculates the sum of the square residuals between two supernova flux 
       arrays usinig a given scale s. Returns the sum.
    """
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
    return np.sum(sq_res)
    
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
    
    ivar_mask[no_data] = 1
    
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
                

def average_new(SN_Array, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask):
    """Modifies the template supernova to be the inverse variance weighted average
       of the scaled data. Returns the new template supernova. 
    """
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
    template.ivar = 1/np.sum(ivars, axis=0)
    template.ivar[no_data] = 0
    template.name = "Composite Spectrum"
    return template
    
#averages with weights of the inverse variances in the spectra
def average(SN_Array, template, medmean):
    """Modifies the template supernova to be the inverse variance weighted average
       of the scaled data. Returns the new template supernova. 
    """
    #create 2D arrays of ell available data
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
    
  
    for i in range(len(dm15_ivars)):
        for j in range(len(dm15_ivars[i])):
            if dm15_ivars[i,j] == 0.0:
                dm15_ivars[i,j] = 1.0
                
#        i_frac = dm15_ivars*fluxes
#        dm15_ivars = i_frac/dm15s
                
#        for i in range(len(dm15s)):
#            if np.all(dm15s[i]) == 0:
#                np.delete(dm15s, i)
#                np.delete(dm15_ivars, i)
    
    
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
    
    #finds and stores the variance data of the template
    template.ivar = 1/np.sum(ivars, axis=0)
    template.ivar[no_data] = 0
    template.name = "Composite Spectrum"
    return template

def boot_mask(SN_Array):
    """Do not need to mask other parameters when bootstrapping"""
    fluxes = []
    ivars  = []
    for SN in SN_Array:
        if len(fluxes) == 0:
            fluxes = np.array([SN.flux])
            ivars  = np.array([SN.ivar])
        else:
            try:
                fluxes = np.append(fluxes, np.array([SN.flux]), axis=0)
                ivars  = np.append(ivars, np.array([SN.ivar]), axis=0)
            except ValueError:
                print "This should never happen!"

    #Adding masks for every parameter for consistency and zero compensation
    flux_mask = np.zeros(len(fluxes[0,:]))
    ivar_mask = np.zeros(len(fluxes[0,:]))
    no_data   = np.where(np.sum(ivars, axis = 0)==0)
    ivar_mask[no_data] = 1
    
    return fluxes, ivars, flux_mask, ivar_mask
    
def boot_avg(SN_Array, template, fluxes, ivars, flux_mask, ivar_mask):
    #Add in flux/ivar mask
    fluxes = np.append(fluxes, np.array([flux_mask]), axis=0)
    ivars  = np.append(ivars, np.array([ivar_mask]), axis=0)
    template.flux  = np.average(fluxes, weights=ivars, axis=0)
    return template

def create_composite(SN_Array, template, scales, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                  reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask):  
    """Scales the current SN_Array to the current template creates a new 
       template from the weighted average of the scaled array. Returns the SN_Array, 
       the template, and the scales used.
    """
    SN_Array, scales = optimize_scales(SN_Array, template)
    template = average_new(SN_Array, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                           reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
        
    return SN_Array, template, scales
    
def do_things(SN_Array, template, scales, medmean, boot, fluxes, ivars, flux_mask, ivar_mask):
    SN_Array, scales = optimize_scales(SN_Array, template)
    if boot:
        template = boot_avg(SN_Array, template, fluxes, ivars, flux_mask, ivar_mask)
    else:
        template = average(SN_Array, template, medmean)
        
    return SN_Array, template, scales
        
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
    
        
#    for i in range(len(strap_matrix)):
#        boot_arr.append([])
#        for j in range(len(strap_matrix[i])):
#            boot_arr[i].append(copy.copy(SN_Array[strap_matrix[i,j]]))
    
    for i in range(len(strap_matrix)):
        boot_arr.append([])
        for j in range(len(strap_matrix[i])):
            boot_arr[i].append(SN_Array[strap_matrix[i,j]])
            

    for p in range(len(boot_arr)):
        lengths = []
        for SN in boot_arr[p]:
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        boot_temp = [SN for SN in boot_arr[p] if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        boot_temp = copy.copy(boot_temp[0])
        
#        fluxes, ivars, flux_mask, ivar_mask = boot_mask(boot_arr[p])
        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
        flux_mask, ivar_mask, dm15_mask, red_mask) = mask(boot_arr[p], boot)
        
        for x in range(iters):
#            SN_Array, template, scales = do_things(boot_arr[p], boot_temp, scales, 1, boot, fluxes, ivars, flux_mask, ivar_mask)
            SN_Array, template, scales = create_composite(boot_arr[p], boot_temp, scales, 1, 
                                                       boot, fluxes, ivars, dm15_ivars, 
                                                       red_ivars, reds, phases, ages, 
                                                       vels, dm15s, flux_mask, ivar_mask, 
                                                       dm15_mask, red_mask)
        boots.append(copy.copy(template))
        
        for i in range(len(SN_Array)):
            SN_Array[i].flux = cpy_array[i].flux
            SN_Array[i].ivar = cpy_array[i].ivar

    
    print "scaling boots..."
    optimize_scales(boots, og_template)
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
    low_ind = np.round(samples * low_pc).astype(int)
    up_ind = np.round(samples * up_pc).astype(int)

    for SN in boots:
        resid.append(SN.flux - og_template.flux)
        
    resid_trans = np.transpose(resid)
    resid_sort = np.sort(resid_trans)
    arr = []
    for elem in resid_sort:
        for i in range (len(elem)):
            arr.append(elem[i])

    plt.hist(arr, 100)
    plt.show()
    
    low_arr = []
    up_arr = []
    for i in range(len(resid_sort)):
        low_arr.append(og_template.flux[i] + resid_sort[i][low_ind])
        up_arr.append(og_template.flux[i] + resid_sort[i][up_ind])
    
    plt.plot(og_template.wavelength, og_template.flux, 'k', linewidth = 4)
    plt.fill_between(og_template.wavelength, low_arr, up_arr, color = 'green')
    plt.show()
    
    return low_arr, up_arr        
    

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
    for SN in SN_Array:
        print SN.name, SN.filename
        SN.flux = SN.flux*1e14
        SN.ivar /= (1e14)**2
        lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
    temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
    try:
        composite = temp[0]
    except IndexError:
        print "No spectra found"
        exit()
    
    spec_bin = spectra_per_bin(SN_Array)

    #finds range of useable data
    template = supernova()
    template = copy.copy(composite)

    #creates main composite
    i = 0
    scales  = []
    iters = 10
    iters_comp = 10
    boot = False
    
    bootstrap = raw_input("Bootstrap? (y/n): ")
    print "Creating composite..."
    
    for i in range(iters_comp):
        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
         flux_mask, ivar_mask, dm15_mask, red_mask) = mask(SN_Array, boot)
#            SN_Array, template, scales = do_things(SN_Array, template, scales, medmean, boot, fluxes, ivars, flux_mask, ivar_mask)
        SN_Array, template, scales = create_composite(SN_Array, template, scales, 
                                                   medmean, boot,fluxes, ivars, 
                                                   dm15_ivars, red_ivars, reds, 
                                                   phases, ages, vels, dm15s, 
                                                   flux_mask, ivar_mask, dm15_mask, 
                                                   red_mask)
    print "Done."
    
    #plot composite with the scaled spectra
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
