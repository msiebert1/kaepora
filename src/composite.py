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
# con = sq3.connect('../data/SNe_3.db')
con = sq3.connect('../data/SNe_14.db')
cur = con.cursor()


def grab(sql_input):
    """Pulls in all columns from the database for the selected query. 
       Replaces all NaN values with 0. Returns the array of supernova objects 
       with the newly added attributes.
    """
    print "Collecting data..."
    SN_Array = []
    # multi_epoch = raw_input("Include multiple epochs? (y/n): ")
    # if multi_epoch == 'y':
    #     multi_epoch = True
    # else:
    #     multi_epoch = False

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
        # phot         = msg.unpackb(row[17])
        phot         = row[17]
        SN.phot      = phot
        SN.low_conf  = []
        SN.up_conf   = []
        SN.spec_bin  = []

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

            # min_phase = events[0]
            # for e in events:
            #     if abs(e.phase) < abs(min_phase.phase):
            #         min_phase = e
            # new_SN_Array.append(min_phase)

            # min_snr = events[0]
            # for e in events:
            #     if abs(e.SNR) > abs(min_snr.SNR):
            #         min_snr = e
            # new_SN_Array.append(min_snr)

            max_range = events[0]
            for e in events:
                if (e.maxwave - e.minwave) > (max_range.maxwave - max_range.minwave):
                    max_range = e
            new_SN_Array.append(max_range)

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
            SN.x2 += 1
        else:
            SN.x1 = 0
            SN.x2 = 0

        # SN.ivar[SN.x1:SN.x1 + 25] = 0.
        # SN.ivar[SN.x2 - 25:SN.x2 + 1] = 0.
        SN.x1 = SN.x1 + 25
        SN.x2 = SN.x2 - 25

                    
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
        # guess = 1.0
        guess = np.average(template.flux[template.x1:template.x2])/np.average(uSN.flux[uSN.x1:uSN.x2])
        # guess = template.flux[2000]/uSN.flux[2000] #this is temporary it does not work for every spectrum 
        if uSN.filename != template.filename:
            u = opt.minimize(sq_residuals, guess, args = (uSN, template, initial), 
                             method = 'Nelder-Mead').x
            scales.append(u)
        else:
            scales.append(1.0)
        
    for i in range(len(unique_arr)):
        unique_arr[i].flux = scales[i]*unique_arr[i].flux
        unique_arr[i].ivar /= (scales[i])**2
        if len(unique_arr[i].low_conf) > 0 and len(unique_arr[i].up_conf) > 0:
            unique_arr[i].low_conf = scales[i]*unique_arr[i].low_conf
            unique_arr[i].up_conf = scales[i]*unique_arr[i].up_conf


        
    return SN_Array, scales
    
    
def sq_residuals(s,SN,comp, initial):
    """Calculates the sum of the square residuals between two supernova flux 
       arrays usinig a given scale s. Returns the sum.
    """
#    print s
    s = np.absolute(s)
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
        # print "no overlap"
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
    temp_fluxes = []
    for i in range(len(SN_Array)):
        fluxes[i] = SN_Array[i].flux
        ivars[i] = SN_Array[i].ivar
        if medmean == 2:
            fluxes[i][fluxes[i] == 0] = np.nan
            if not boot:
                ages[i][ages[i] == 0] = np.nan
                vels[i][vels[i] == 0] = np.nan
                dm15s[i][dm15s[i] == 0] = np.nan
                reds[i][reds[i] == 0] = np.nan
        # temp_fluxes.append(np.ma.masked_where(fluxes[i] == 0, fluxes[i]))    
        # temp_fluxes.append(fluxes[i][fluxes[i] == 0] = np.nan) 
    
    # print temp_fluxes

    if medmean == 1: 
        template.flux  = np.average(fluxes, weights=ivars, axis=0)
        if not boot:
            template.phase_array   = np.average(ages, weights=ivars, axis=0)
            template.vel   = np.average(vels, weights=ivars, axis=0)
            template.dm15  = np.average(dm15s, weights=dm15_ivars, axis=0)
            template.red_array = np.average(np.array(reds), weights = red_ivars, axis=0)
    if medmean == 2:
        # template.flux  = np.median(fluxes, axis=0)
        # template.flux  = np.ma.median(temp_fluxes, axis=0).filled(0)
        template.flux = np.nanmedian(fluxes, axis=0)
        if not boot:
            template.phase_array   = np.nanmedian(ages, axis=0)
            template.vel   = np.nanmedian(vels, axis=0)
            template.dm15  = np.nanmedian(dm15s, axis=0)
            template.red_array = np.nanmedian(reds, axis=0)
    
    #finds and stores the variance data of the template
    no_data   = np.where(np.sum(ivars, axis = 0)==0)
    template.ivar = np.sum(ivars, axis=0)
    template.ivar[no_data] = 0
    template.name = "Composite Spectrum"
    return template

        
def bootstrapping (SN_Array, samples, scales, og_template, iters, medmean):
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
        # print p
        lengths = []
        for SN in boot_arr[p]:
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        boot_temp = [SN for SN in boot_arr[p] if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        boot_temp = copy.copy(boot_temp[0])
        
        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
         flux_mask, ivar_mask, dm15_mask, red_mask) = mask(boot_arr[p], boot)
        for x in range(iters):
            SN_Array, scales = optimize_scales(boot_arr[p], boot_temp, False)
            template = average(boot_arr[p], boot_temp, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
                               reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
        boots.append(copy.copy(template))

    print "scaling boots..."
    temp1, scales = optimize_scales(boots, og_template, True)
    if medmean == 1:
        optimize_scales(boots, og_template, False)

    # print "plotting..."
    # for SN in boots:
    #     plt.plot(SN.wavelength, SN.flux, 'g')
    # plt.plot(og_template.wavelength,og_template.flux, 'k', linewidth = 4)
    # plt.show()
    
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

    # plt.hist(arr, 100)
    # plt.show()
    
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
    
    # plt.plot(og_template.wavelength, og_template.flux, 'k', linewidth = 4)
    # plt.fill_between(og_template.wavelength, low_arr, up_arr, color = 'green')
    # plt.show()
    
    return np.asarray(low_arr), np.asarray(up_arr)        
    
def is_bad_data(SN, bad_files, bad_ivars):
    for el in bad_files:
        if SN.filename == el:
            print 'bad file removed'
            return True
    for el in bad_ivars:
        if SN.filename == el:
            return True
    return False

def remove_peculiars(SN_Array, file):
    SN_Array_no_pecs = []
    with open(file) as f:
        names = np.loadtxt(f, dtype = str)
        for SN in SN_Array:
            if SN.name not in names:
                SN_Array_no_pecs.append(SN)

    return SN_Array_no_pecs

def build_av_dict(file):
     with open(file) as f:
        lines = f.readlines()

        av_dict = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                av_dict[l[1].lower()] = float(l[18])

     return av_dict

def split_list(n):
    """will return the list index"""
    return [(x+1) for x,y in zip(n, n[1:]) if y-x != 1.]

def get_sub_list(my_list):
    """will split the list base on the index"""
    my_index = split_list(my_list)
    output = list()
    prev = 0
    for index in my_index:
        new_list = [ x for x in my_list[prev:] if x < index]
        output.append(new_list)
        prev += len(new_list)
    output.append([ x for x in my_list[prev:]])
    return output

def fix_weights(SN_Array):

    for SN in SN_Array:
        zero_inds = np.where(SN.ivar == 0)[0]
        zero_ranges = get_sub_list(zero_inds)

        for r in zero_ranges:
            if r[0] > SN.x1 and r[-1] < SN.x2:
                SN.ivar[r] = (SN.ivar[r[0]-1] + SN.ivar[r[-1]+1])/2.
                SN.flux[r] = (SN.flux[r[0]-1] + SN.flux[r[-1]+1])/2.
    
    
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
    bad_files =        ['SN05hc_051018_r01_NTT_EM.dat', '2003du_20030501_4066_11015_00.dat',
                        '2002er_20020901_3213_9175_00.dat', '2003du_20030503_4070_10994_00.dat',
                        '2003du_20030512_4070_11015_00.dat', 'sn2006nz-20061124-ntt.dat',
                        'sn2001eh-20010925.757-hst.flm', 'sn2001ep-20011102.871-hst.flm',
                        '2005cf_20050601_3243_9720_00.dat', 'sn2006cm-20060529.46-fast.flm',
                        '2002bo_20020321_3357_7725_00.dat', 'sn1994S-19940612.26-mmt.flm',
                        'sn2003cg-20030329.27-mmt.flm', 'sn1995ac-19950929.27-fast.flm',
                        'sn2007hj-20070903.28-fast.flm', '2000E_20000131_3274_7349_00.dat',
                        'sn2006cj-20060521.29-fast.flm', 'sn2006oa-20061125.08-fast.flm',
                        'sn2005cf-20050609.5-uvot-clip.flm', 'sn2006kf-20061030.385-ui.flm',
                        'SN07bd_070412_b01_DUP_WF.dat', 'SN09ad_090223_b01_DUP_WF.dat',
                        'SN05kc_051124_b01_DUP_MS.dat', 'sn2005eq-20051011.376-ui-corrected.flm',
                        'sn2006et-20060919.345-ui.flm', 'sn2007cq-20070623.431-ui.flm',
                        'sn1997bq-19970408.14-mmt.flm', 'sn2006lf-20061028.51-fast.flm',
                        'sn2005eq-20051002.51-fast.flm', 'sn1995bd-19951223.34-fast.flm',
                        'sn1998ab-19980403.38-fast.flm', 'sn1994M-19940612.22-mmt.flm',
                        '2006X_20060209_3834_8139_00.dat', '2003du_20030508_4066_10997_00.dat',
                        'SN05ke_051125_b01_T60_CS.dat', 'sn1994s-19940616-uoi.flm',
                        '1994D_19940317_2999_10549_00.dat', '2005cf_20050608_3365_9997_00.dat',
                        '2002dj_20020620_3210_9231_00.dat', 'sn2006lf-20061029.40-fast.flm',
                        'sn2006lf-20061030.41-fast.flm', 'sn2006lf-20061031.42-fast.flm',
                        'sn2006lf-20061101.40-fast.flm', 'sn2006lf-20061111.36-fast.flm',
                        'sn2006lf-20061112.37-fast.flm','sn2006lf-20061113.36-fast.flm',
                        'sn2006lf-20061115.42-fast.flm','sn2006lf-20061116.39-fast.flm',
                        'sn2006lf-20061117.42-fast.flm','sn2006lf-20061119.35-fast.flm',
                        'sn2006lf-20061122.35-fast.flm','sn2006lf-20061125.34-fast.flm',
                        'sn2006lf-20061214.28-fast.flm','sn2006lf-20061216.30-fast.flm',
                        'sn2006lf-20061226.26-fast.flm','sn2006lf-20061227.29-fast.flm',
                        'sn1995al-19951114.52-fast.flm', '2002er_20020901_3213_9175_00.dat',
                        'sn2003kf-20031216.37-fast.flm','SN05bg_050419_b01_DUP_WF.dat',
                        '2003du_20040202_3182_9333_00.dat', '2002bo_20020328_3099_8801_00.dat',
                        '2005cf_20050603_3721_8786_00.dat', 'SN06hx_061005_g01_NTT_EM.dat',
                        '2003du_20030429_3428_9436_00.dat', '2000E_20000127_3213_7513_00.dat',
                        'sn2006bt-20060427.629-br.flm','sn2003ic-20030919.40-fast.flm',
                        '2000E_20000130_3274_7356_00.dat','sn2006cm-20060528.45-fast.flm',
                        'sn2006nz-20061117.18-fast.flm','sn1996ai-19960620-uo.flm',
                        'sn2006cm-20060528.424-ui.flm','2003du_20030506_4068_10991_00.dat',
                        'sn2005A-20050111.04-ldss2.flm','sn2003iv-20031023.469-ui.flm',
                        'SN07al_070314_b01_DUP_BC.dat', 'sn2006ke-20061024.397-ui.flm',
                        'sn2002bf-20020307-os.flm','sn1996bk-19961015.10-fast.flm',
                        'SN06mr_061113_r01_BAA_IM.dat','sn2003ic-20030927.38-fast.flm',
                        'SN08hv_081220_b01_DUP_WF.dat','sn2002de-20020614.19-fast.flm',
                        'sn2003ic-20030929.35-fast.flm','sn2006cj-20060529.26-fast.flm',
                        'SN06mr_061116_b01_DUP_WF.dat','SN07al_070319_b01_DUP_BC.dat',
                        '2003du_20030515_4969_9241_00.dat','SN07N_070131_b01_NTT_EM.dat',
                        'sn2006ke-20061030.525-ui.flm','sn2003ic-20031002.25-fast.flm',
                        'sn2005hf-20051027.27-fast.flm', 'SN06mr_061119_b01_DUP_WF.dat',
                        'sn2005hf-20051028.22-fast.flm','SN06mr_061119_b01_HIL_BC.dat',
                        'sn2001ic-20011211-ui.flm','sn2005mz-20060121.13-fast.flm',
                        'sn2000cn-20000623.34-fast.flm','sn2007al-20070320.11-ldss3.flm',
                        'sn2007ax-20070407.16-fast.flm','SN07on_071125_b01_CLA_MA.dat', 
                        'SN08hv_081227_b01_DUP_WF.dat','sn2005hf-20051030.22-fast.flm',
                        'sn2007al-20070321.33-fast.flm','sn1999gh-19991212.50-fast.flm',
                        'sn2001N-20010202.47-fast.flm','sn2002bf-20020316.32-fast.flm',
                        'SN06mr_061122_b01_DUP_WF.dat','sn2005hf-20051031.24-fast.flm',
                        'sn1998bp-19980516.44-fast.flm','sn2003hu-20031002.13-fast.flm',
                        'sn2005hf-20051101.22-fast.flm','SN06bd_060330_b01_DUP_WF.dat',
                        'sn2007ci-20070609.244-ui-corrected.flm','SN07al_070326_b01_DUP_BC.dat',
                        'sn1998bp-19980518.36-fast.flm','SN05kc_051208_b01_T60_CS.dat',
                        'sn1999gh-19991217-ui.flm','sn2007bc-20070426.387-ui-corrected.flm',
                        'sn2002do-20020704.39-fast.flm','SN07on_071203_b01_DUP_WF.dat',
                        'sn2005mz-20060129.11-fast.flm','SN06gt_061013_b01_DUP_WF.dat',
                        'sn1998bp-19980522-r.flm','sn2005hf-20051106.34-fast.flm',
                        'sn2006cm-20060620.410-ui.flm','sn2000cn-20000704.32-fast.flm',
                        'SN08R_080219feb08_b01_CLA_MA.dat','SN07jg_071016_b01_DUP_BC.dat',
                        '2003du_20030530_4060_10974_00.dat','sn2006nz-20061213-ntt.dat',
                        'sn2006te-20070126.326-ui.flm','SN08ia_090122_b01_CLA_LD.dat',
                        'sn2001G-20010423.19-fast.flm','SN06X_060524_b01_DUP_BC.dat',
                        'sn2006X-20060221.40-fast.flm','2006X_20060219_3731_8515_00.dat',
                        '2006X_20060221_3981_8865_00.dat','sn1999cl-19990614.18-fast.flm',
                        'sn2006X-20060222.41-fast.flm','sn2006x-20060222.413-ui.flm',
                        '2006X_20060225_3734_8223_00.dat','sn2006X-20060225.36-fast.flm',
                        '2006X_20060227_3918_8203_00.dat','sn2006X-20060227.44-fast.flm',
                        'sn2006X-20060228.34-fast.flm','2002er_20020916_3336_8734_00.dat',
                        'sn2006X-20060302.47-fast.flm','2006X_20060304_3783_8272_00.dat',
                        'sn2006X-20060304.51-fast.flm','2006X_20060307_3861_8130_00.dat',
                        'sn2006X-20060309.30-fast.flm','2002er_20020926_3489_8768_00.dat',
                        '2002er_20021010_3560_9363_00.dat','1996X_19960613_2806_10203_00.dat']
                        #SN05hc_051018_r01_NTT_EM.dat very noisy
                        #2003du_20030501_4066_11015_00.dat very large negative value
                        #2002er_20020901_3213_9175_00.dat gap in spectrum
                        #2003du_20030503_4070_10994_00.dat very strong emission line?
                        #2003du_20030512_4070_11015_00.dat very negative values
                        #sn2006nz-20061124-ntt.dat very negative values
                        #sn2001eh-20010925.757-hst.flm some interpolated sections
                        #sn2001ep-20011102.871-hst.flm causes problem not sure why
                        #2005cf_20050601_3243_9720_00.dat causes problem not sure why
                        #sn2006cm-20060529.46-fast.flm weird snr
                        #2002bo_20020321_3357_7725_00.dat doesn't scale properly, flux very small at low wavelength
                        #sn1994S-19940612.26-mmt.flm silicon line was interpolated
                        #sn2003cg-20030329.27-mmt.flm SNR above 600 biases averaging
                        #sn1995ac-19950929.27-fast.flm noisy
                        #sn2007hj-20070903.28-fast.flm some interpolated sections
                        #2000E_20000131_3274_7349_00.dat spectrum seems to be blueshifted
                        #sn2006cj-20060521.29-fast.flm very noisy and some interpolation
                        #sn2006oa-20061125.08-fast.flm some interpolated sections
                        #sn2005cf-20050609.5-uvot-clip.flm high snr drags composite down
                        #sn2006kf-20061030.385-ui.flm some interpolated sections
                        #SN07bd_070412_b01_DUP_WF.dat some interpolated sections
                        #SN09ad_090223_b01_DUP_WF.dat some interpolated sections
                        #SN05kc_051124_b01_DUP_MS.dat some interpolated sections
                        #sn2006et-20060919.345-ui.flm some interpolated sections
                        #sn2007cq-20070623.431-ui.flm some interpolated sections
                        #sn1997bq-19970408.14-mmt.flm telluric?
                        #sn2006lf-20061028.51-fast.flm some interpolated sections, variance spectrum seems wrong
                        #sn2005eq-20051002.51-fast.flm some interpolated sections
                        #sn1995bd-19951223.34-fast.flm some interpolated sections
                        #sn1998ab-19980403.38-fast.flm some interpolated sections
                        #sn1994M-19940612.22-mmt.flm large interpolated section
                        #2006X_20060209_3834_8139_00.dat host correction causes
                        #2003du_20030508_4066_10997_00.dat very large negative value
                        #SN05ke_051125_b01_T60_CS.dat some interpolated sections
                        #sn1994s-19940616-uoi.flm some interpolated sections
                        #1994D_19940317_2999_10549_00.dat not joined properly, telluric absorption
                        #2005cf_20050608_3365_9997_00.dat telluric absorption
                        #2002dj_20020620_3210_9231_00.dat telluric absorption
                        #2002er_20020901_3213_9175_00.dat large gap
                        #
                        #All 2006lf cfa data seems to have large ivar (skews data)
                        #sn2003kf-20031216.37-fast.flm seems to have large ivar
                        #SN05bg_050419_b01_DUP_WF.dat ivar very large at large wavelengths
                        #2003du_20040202_3182_9333_00.dat variance blows up
                        #2002bo_20020328_3099_8801_00.dat large interpolated section
                        #2005cf_20050603_3721_8786_00.dat large interpolated section
                        #SN06hx_061005_g01_NTT_EM.dat slope is off after MW correction
                        #2003du_20030429_3428_9436_00.dat telluric absorption
                        #2000E_20000127_3213_7513_00.dat weird slope and noise
                        #sn2006bt-20060427.629-br.flm slope is off after MW correction
                        #sn2003ic-20030919.40-fast.flm slope is off after MW correction
                        #2000E_20000130_3274_7356_00.dat slope is off after MW correction
                        #sn2006cm-20060528.45-fast.flm slope is off after MW correction
                        #sn2006nz-20061117.18-fast.flm slope is off after MW correction
                        #sn1996ai-19960620-uo.flm large interpolated section
                        #sn2006cm-20060528.424-ui.flm slope is off affter MW correction
                        #2003du_20030506_4068_10991_00.dat large negative value
                        #sn2005A-20050111.04-ldss2.flm slope is off after MW correction
                        #sn2003iv-20031023.469-ui.flm slope is off after MW correction
                        #SN07al_070314_b01_DUP_BC.dat slope is off after MW correction
                        #sn2006ke-20061024.397-ui.flm slope is off after MW correction
                        #sn2002bf-20020307-os.flm large blueshift?
                        #sn1996bk-19961015.10-fast.flm slope is off after MW correction
                        #SN06mr_061113_r01_BAA_IM.dat weird shape
                        #sn2003ic-20030927.38-fast.flm slope is off after MW correction
                        #SN08hv_081220_b01_DUP_WF.dat slope is off after MW correction
                        #sn2002de-20020614.19-fast.flm slope is off after MW correction
                        #sn2003ic-20030929.35-fast.flm slope is off after MW correction
                        #sn2006cj-20060529.26-fast.flm slope is off after MW correction
                        #SN06mr_061116_b01_DUP_WF.dat continuum seems off
                        #SN07al_070319_b01_DUP_BC.dat continuum seems off
                        #2003du_20030515_4969_9241_00.dat large blueshift?
                        #SN07N_070131_b01_NTT_EM.dat continuum seems off
                        #sn2006ke-20061030.525-ui.flm continuum seems off
                        #sn2003ic-20031002.25-fast.flm slope is off after MW correction
                        #sn2005hf-20051027.27-fast.flm slope is off after MW correction
                        #SN06mr_061119_b01_DUP_WF.dat continuum seems off
                        #sn2005hf-20051028.22-fast.flm slope is off after MW correction
                        #SN06mr_061119_b01_HIL_BC.dat continuum seems off
                        #sn2001ic-20011211-ui.flm almost entirely flat
                        #sn2005mz-20060121.13-fast.flm slope is off after MW correction
                        #sn2000cn-20000623.34-fast.flm slope is off after MW correction
                        #sn2007al-20070320.11-ldss3.flm continuum seems off
                        #sn2007ax-20070407.16-fast.flm slope is off after MW correction
                        #SN07on_071125_b01_CLA_MA.dat continuum seems off, telluric contamination
                        #SN08hv_081227_b01_DUP_WF.dat continuum seems off
                        #sn2005hf-20051030.22-fast.flm spectrum is flat
                        #sn2007al-20070321.33-fast.flm continuum seems off
                        #sn1999gh-19991212.50-fast.flm continuum seems off
                        #sn2001N-20010202.47-fast.flm slope is off of MW correction
                        #sn2002bf-20020316.32-fast.flm slope is off of MW correction
                        #SN06mr_061122_b01_DUP_WF.dat continuum seems off
                        #sn2005hf-20051031.24-fast.flm spectrum is flat
                        #sn1998bp-19980516.44-fast.flm continuum seems off
                        #sn2003hu-20031002.13-fast.flm spectrum is flat
                        #sn2005hf-20051101.22-fast.flm continuum seems off
                        #SN06bd_060330_b01_DUP_WF.dat continuum seems off
                        #sn2007ci-20070609.244-ui-corrected.flm continuum seems off
                        #SN07al_070326_b01_DUP_BC.dat continuum seems off
                        #sn1998bp-19980518.36-fast.flm continuum seems off
                        #SN05kc_051208_b01_T60_CS.dat variance spectrum is wrong
                        #sn1999gh-19991217-ui.flm continuum seems off 
                        #sn2007bc-20070426.387-ui-corrected.flm continuum seems off
                        #sn2002do-20020704.39-fast.flm spectrum is flat
                        #SN07on_071203_b01_DUP_WF.dat continuum seems off
                        #sn2005mz-20060129.11-fast.flm continuum seems off
                        #SN06gt_061013_b01_DUP_WF.dat continuum seems off
                        #sn1998bp-19980522-r.flm continuum seems off
                        #sn2005hf-20051106.34-fast.flm spectrum is flat
                        #sn2006cm-20060620.410-ui.flm continuum seems off
                        #sn2000cn-20000704.32-fast.flm continuum seems off
                        #SN08R_080219feb08_b01_CLA_MA.dat continuum seems off, telluric absorption
                        #SN07jg_071016_b01_DUP_BC.dat very noisy, many interpolated sections
                        #2003du_20030530_4060_10974_00.dat large negative values
                        #sn2006nz-20061213-ntt.dat large negative values
                        #sn2006te-20070126.326-ui.flm spectrum is flat
                        #SN08ia_090122_b01_CLA_LD.dat redshifted?
                        #sn2001G-20010423.19-fast.flm continuum seems off
                        #SN06X_060524_b01_DUP_BC.dat doesnt scale properly
                        #sn2006X-20060221.40-fast.flm blueshifted?
                        #2006X_20060219_3731_8515_00.dat blueshifted?
                        #2006X_20060221_3981_8865_00.dat blueshifted?
                        #sn1999cl-19990614.18-fast.flm variance seems wrong
                        #sn2006X-20060222.41-fast.flm blueshifted?
                        #sn2006x-20060222.413-ui.flm blueshifted?
                        #2006X_20060225_3734_8223_00.dat blueshifted?
                        #sn2006X-20060225.36-fast.flm blueshifted?
                        #2006X_20060227_3918_8203_00.dat blueshifted?
                        #sn2006X-20060227.44-fast.flm blueshifted?
                        #sn2006X-20060228.34-fast.flm blueshifted?
                        #2002er_20020916_3336_8734_00.dat seems off near red edge
                        #sn2006X-20060302.47-fast.flm blueshifted?
                        #2006X_20060304_3783_8272_00.dat blueshifted?
                        #sn2006X-20060304.51-fast.flm blueshifted?
                        #2006X_20060307_3861_8130_00.dat blueshifted?
                        #sn2006X-20060309.30-fast.flm blueshifted?
                        #2002er_20020926_3489_8768_00.dat continuum seems off
                        #2002er_20021010_3560_9363_00.dat continuum seems off
                        #1996X_19960613_2806_10203_00.dat blueshifted?

                        # sn1994d-19940603.flm ??

                        #sn2001az-20010430-ui.flm

                        #Weird continuum
                        #sn1997E-19970116.25-fast.flm   0.1221
                        #sn1996bl-19961018.18-fast.flm  0.2945
                        #sn1998ef-19981024.30-fast.flm  -0.0303
                        #sn2003it-20031019.18-fast.flm  0.0932

    reddened_spectra = []
    
    bad_ivars = []
    for SN in SN_Array:
        if True in np.isnan(SN.ivar):
            bad_ivars.append(SN.filename)

    good_SN_Array = [SN for SN in SN_Array if not is_bad_data(SN, bad_files, bad_ivars)]
    SN_Array = good_SN_Array

    print len(bad_ivars), 'spectra with nan ivars removed'

    #remove peculiar Ias
    SN_Array = remove_peculiars(SN_Array,'../data/info_files/pec_Ias.txt')
    print 'Peculiar Ias removed', len(SN_Array), 'spectra left'

    
    # mags = np.sort(mg.ab_mags(SN_Array), axis = 0)
    # for i in range(len(mags)):
    #     print mags[i][0], mags[i][1] 

    corrected_SNs = []
    av_dict = build_av_dict('../data/info_files/lowz_rv25_all.fitres')
    r_v = 2.5
    for SN in SN_Array:
        # if SN.wavelength[SN.x2] > 10240 and SN.wavelength[SN.x2] < 10270:
        #     print SN.filename, '                                  here'
        print SN.name, SN.filename, SN.phase
        # if np.average(SN.flux[SN.x1:SN.x2]) > 1.e-13:
            # SN.flux = 1.e-15*SN.flux
        pre_scale = (1.e-15/np.average(SN.flux[SN.x1:SN.x2]))
        SN.flux = pre_scale*SN.flux
        SN.ivar = SN.ivar/(pre_scale*pre_scale)
        old_wave = SN.wavelength*u.Angstrom        # wavelengths
        old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
        spec1d = Spectrum1D.from_array(old_wave, old_flux)
        old_ivar = SN.ivar*u.Unit('W m-2 angstrom-1 sr-1')
        # new_flux = test_dered.dered(sne, SN.name, spec1d.wavelength, spec1d.flux)
        new_flux, new_ivar, corrected = test_dered.host_correction(av_dict, r_v, SN.name, old_wave, old_flux, old_ivar)
        SN.flux = new_flux.value
        SN.ivar = new_ivar.value
        # lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        if corrected:
            corrected_SNs.append(SN)
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

    SN_Array = corrected_SNs
    print len(SN_Array), 'SNs with host corrections'


    # fix_weights(SN_Array)


    # good_SNs = []
    # lengths = []
    # for SN in SN_Array:
    #     plt.figure(num = 2, dpi = 100, figsize = [35, 17], facecolor = 'w')
    #     plt.subplot(2,1,1)
    #     plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
    #     plt.subplot(2,1,2)
    #     plt.plot(SN.wavelength[SN.x1:SN.x2], 1./(SN.ivar[SN.x1:SN.x2])**.5)
    #     # plt.plot(SN.wavelength[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2])
    #     plt.show()
    #     print SN.SNR, SN.filename, SN.source
    #     accept = raw_input("Accept? (y/n): ")
    #     if accept is 'y':
    #         print "Added"
    #         good_SNs.append(SN)
    #         lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

    # SN_Array = good_SNs
    # print len(SN_Array)

    temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
    try:
        composite = temp[0]
    except IndexError:
        print "No spectra found"
        exit()

    # spec_bin = spectra_per_bin(SN_Array)

    #finds range of useable data
    template = supernova()
    template = copy.copy(composite)
    template.spec_bin = spectra_per_bin(SN_Array)

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
    # for i in range(len(SN_Array)):
    #     plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
    # plt.plot(template.wavelength, template.flux, 'k', linewidth = 4)
    # plt.show()
    
    #for updating one spectrum at a time
    # num_plots = len(SN_Array)
    # for i in range(num_plots):
    #     sub_sns = copy.copy(SN_Array[0:i+1])
    #     print SN_Array[i].filename, SN_Array[i].phase, SN_Array[i].source, SN_Array[i].SNR
    #     optimize_scales(sub_sns, template, True)
    #     (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
    #      flux_mask, ivar_mask, dm15_mask, red_mask) = mask(sub_sns, boot)
    #     for j in range(iters_comp):
    #         SN_Array, scales = optimize_scales(SN_Array, template, False)
    #         template = average(sub_sns, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
    #                            reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)

    #     plt.figure(num = 2, dpi = 100, figsize = [30, 15], facecolor = 'w')
    #     plt.subplot(2,1,1)
    #     # plt.plot(sub_sns[-1].wavelength, sub_sns[-1].flux)
    #     for SN in sub_sns:
    #         plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
    #     plt.plot(template.wavelength[SN.x1:SN.x2], template.flux[SN.x1:SN.x2], 'k', linewidth = 4)
    #     plt.subplot(2,1,2)
    #     # plt.plot(sub_sns[-1].wavelength, sub_sns[-1].ivar)
    #     for SN in sub_sns:
    #         plt.plot(SN.wavelength[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2])
    #     plt.show()


    # bootstrap = raw_input("Bootstrap? (y/n): ")
    bootstrap = 'y'
    # bootstrap = 'y'
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
        # pass
        for i in range(len(SN_Array)):
            plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
        plt.plot(template.wavelength, template.flux, 'k', linewidth = 4)
        plt.show()

    #create bootstrap composites
    else:
        scales  = []
        print "Bootstrapping"
        # samples = int (raw_input("# of samples:"))
        samples = 100
        # low_conf, up_conf = bootstrapping(SN_Array, samples, scales, template, iters, medmean)
        template.low_conf, template.up_conf = bootstrapping(SN_Array, samples, scales, template, iters, medmean)
    
    non_zero_data = np.where(template.flux != 0)
    non_zero_data = np.array(non_zero_data[0])
    if len(non_zero_data) > 0:
        SN.x1 = non_zero_data[0]
        SN.x2 = non_zero_data[-1]
        SN.x2 += 1

    template.ivar = 1./template.ivar
    template.ivar[0:template.x1] = 0.
    template.ivar[template.x2:] = 0.

    # table = Table([template.wavelength, template.flux, template.ivar, template.phase_array, 
    #                template.vel, template.dm15, template.red_array, low_conf, up_conf, spec_bin], 
    #                names = ('Wavelength', 'Flux', 'Variance', 'Age', 'Velocity', 
    #                         'Dm_15', 'Redshift', 'Lower Confidence', 'Upper Confidence', 'Spectra Per Bin'))
    # if save_file == 'y':
    #     table.write(template.savedname, format='ascii')

    # return table
    return template

if __name__ == "__main__":
    main()
