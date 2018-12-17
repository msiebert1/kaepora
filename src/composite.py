"""
Spectra composite program
Authors: Sam, Yixian, Aaron
"""
#example: python Run.py nb "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN and phase between -1 and 1"

#version specific dependencies:
#msgpack_python version 0.4.6
#msgpack_numpy version 0.3.6

#python packages
import matplotlib.pyplot as plt
import numpy as np
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
from scipy.special import erf
import scipy.optimize as opt
import copy
from astropy import units as u
from specutils import Spectrum1D

#routines in repository:
import test_dered
import questionable_spectra as qspec
import telluric_spectra as tspec
import query_db as qdb
import spectral_analysis as sa
import gini

np.set_printoptions(threshold=np.nan)
mn.patch()

class supernova(object):
    """A generic class to represent a spectrum and its associated metadata
    """
    def __init__(self, wavelength = None, flux = None, low_conf=[], up_conf=[], ivar = None):
        if wavelength is not None:
            self.wavelength = wavelength
            self.flux = flux
            if ivar is None:
                self.ivar = np.zeros(len(flux))
            else:
                self.ivar = ivar
            self.low_conf = low_conf
            self.up_conf =up_conf
            self.x1 = 0
            self.x2 = len(wavelength) - 1

def store_phot_data(SN, row):
    """Assigns attributes to a supernova object from the Photometry table in the
    SQL database. 

    Args:
        SN: A supernova object that has been initialized by the 'grab' function.
        row: The row from the Photometry table associated with this SN object.

    Yields:
        A SN object with its event specific attributes.

    Note:
        This function is only called if 'join' is in the in the SQL query. This
        indicates that the user wants event specific metadata. 
    """

    phot_row = row[19:]
    SN.ra = phot_row[1]
    SN.dec = phot_row[2]
    SN.zCMB_salt, SN.e_zCMB_salt, SN.Bmag_salt, SN.e_Bmag_salt, SN.s_salt, SN.e_s_salt, \
    SN.c_salt, SN.e_c_salt, SN.mu_salt, SN.e_mu_salt = phot_row[3:13]
    SN.zCMB_salt2, SN.e_zCMB_salt2, SN.Bmag_salt2, SN.e_Bmag_salt2, SN.x1_salt2, \
    SN.e_x1_salt2, SN.c_salt2, SN.e_c_salt2, SN.mu_salt2, SN.e_mu_salt2 = phot_row[13:23]
    SN.zCMB_mlcs31, SN.e_zCMB_mlcs31, SN.mu_mlcs31, SN.e_mu_mlcs31, SN.delta_mlcs31, \
    SN.e_delta_mlcs31, SN.av_mlcs31, SN.e_av_mlcs31 = phot_row[23:31]
    SN.zCMB_mlcs17, SN.e_zCMB_mlcs17, SN.mu_mlcs17, SN.e_mu_mlcs17, SN.delta_mlcs17, \
    SN.e_delta_mlcs17, SN.av_mlcs17, SN.e_av_mlcs17 = phot_row[31:39]
    SN.glon_host, SN.glat_host, SN.cz_host, SN.czLG_host, SN.czCMB_host, SN.mtype_host, \
    SN.xpos_host, SN.ypos_host, SN.t1_host, SN.filt_host, SN.Ebv_host = phot_row[39:50]
    SN.zCMB_lc, SN.zhel_lc, SN.mb_lc, SN.e_mb_lc, SN.c_lc, SN.e_c_lc, SN.x1_lc, \
    SN.e_x1_lc, SN.logMst_lc, SN.e_logMst_lc, SN.tmax_lc, SN.e_tmax_lc, \
    SN.cov_mb_s_lc, SN.cov_mb_c_lc, SN.cov_s_c_lc, SN.bias_lc = phot_row[50:66]
    SN.av_25 = phot_row[66]
    SN.dm15_source = phot_row[67]
    SN.dm15_from_fits = phot_row[68]
    SN.e_dm15 = phot_row[69]
    SN.sep = phot_row[70]
    SN.ned_host = phot_row[71]
    SN.v_at_max = phot_row[72]
    SN.e_v = phot_row[73]
    SN.light_curves = msg.unpackb(phot_row[-2])
    SN.csp_light_curves = msg.unpackb(phot_row[-1])

def grab(sql_input, multi_epoch = False, make_corr = True, 
         selection = 'max_coverage', grab_all=False):
    """A primary function for interacting with the database. The user specifies 
    a desired subsample via an SQL query. Spectra are stored as supernova objects
    with their associated metadata as attributes.

    Args:
        sql_input: The SQL query string 

    Keyword arguments:
        multi_epoch: If True, include all spectra for a given SN that satisify 
            the query. If False, choose one 1 spectrum per SN based on the 
            selection keyword.
        make_corr: If True, remove spectra that have been marked as 
            'questionable', peculiar events (Iax), and spectra that do not have 
            host extinction estimates.
        selection: If multi_epoch is False, this string defines the selection
            criteria for choosing a single spectrum from a SN. Options are:
            'maximum_coverage'(default): largest wavelength range
            'maximum_coverage_choose_uv': largest wavelength range but prioritize
                hst and swift spectra
            'choose_bluest': smallest minimum wavelength
            'max_snr': highest signal to noise
            'accurate_phase': closest to middle of the phase bin. TODO: implement
                this without parsing the query (currently requires uncommenting
                code)
            'max_coverage_splice': allows multiple spectra from the same SN as 
                long as overlap is < 500 A
        grab_all: If True, ignore other arguments and return all data that 
            satisfy the SQL query.

    Returns:
        An array of supernova objects populated with metadata retrieved from 
        the SQL query.

    """
    print "Collecting data..."
    #Connect to database
    #Make sure your database file is in this location
    con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
    cur = con.cursor()

    SN_Array = []

    get_phot = False
    if "join" in sql_input:
        get_phot = True

    cur.execute(sql_input)

    #UNCOMMMENT if you want to parse query for average phase (assumes a syntax)
    # if 'phase' in sql_input:
    #   split_query = sql_input.split()
    #   p_index = split_query.index('phase')
    #   p_low = float(split_query[p_index+2])
    #   p_high = float(split_query[p_index+6])
    #   p_avg = (p_low+p_high)/2.

    for row in cur:
        """retrieve spectral metadata. TODO: Move or delete event specific metadata
        from the "Supernovae" table.
        """

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
        SN.B_minus_V = row[9]
        SN.velocity  = row[10]
        SN.morph     = row[11]
        SN.carbon    = row[12]
        SN.GasRich   = row[13]
        SN.SNR       = row[14]
        SN.resid     = row[15]
        interp       = msg.unpackb(row[16])
        SN.interp    = interp
        SN.mjd         = row[17]
        SN.ref       = row[18]

        #retrieve event specific metadata if desired
        if get_phot:
            store_phot_data(SN, row)

        SN.low_conf  = []
        SN.up_conf   = []
        SN.spec_bin  = []

        try:
            SN.wavelength = SN.interp[0,:]
            SN.flux       = SN.interp[1,:]
            SN.ivar       = SN.interp[2,:]
        except TypeError:
            print "ERROR: ", SN.filename, SN.interp
            continue
        SN_Array.append(SN)

    if grab_all:
        return SN_Array

    if make_corr:
        """remove flagged files from questionable_spectra.py, remove spectra
        with bad ivar data, remove spectra without host extinction esimates, and
        remove spectra marked as 'peculiar' (currently only Iax)
        """
        bad_files = qspec.bad_files()

        bad_ivars = []
        for SN in SN_Array:
            # print SN.filename 
            if len(np.where(np.isnan(SN.ivar))[0] == True) == 5500:
                bad_ivars.append(SN.filename)
                # plt.plot(SN.wavelength,SN.flux)
                # plt.show()
                # plt.plot(SN.wavelength,SN.ivar)
                # plt.show()
        if len(bad_ivars) > 0:
            print "Generate variance failed for: ", bad_ivars

        len_before = len(SN_Array)
        good_SN_Array = [SN for SN in SN_Array 
                            if not is_bad_data(SN, bad_files, bad_ivars)]
        SN_Array = good_SN_Array
        print (len_before - len(SN_Array), 'questionable spectra removed', 
               len(SN_Array), 'spectra left')

        # remove peculiar Ias
        len_before = len(SN_Array)
        SN_Array = remove_peculiars(SN_Array,'../data/info_files/pec_Ias.txt')
        print (len_before - len(SN_Array), 'Peculiar Ias removed', 
               len(SN_Array), 'spectra left')

        SN_Array = check_host_corrections(SN_Array)

    if not multi_epoch:
        #choose one spectrum per object based on selection criteria
        bad_files = qspec.bad_files()
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
            if selection == 'max_coverage':
                max_range = events[0]
                for e in events:
                    if (e.maxwave - e.minwave) > (max_range.maxwave - max_range.minwave):
                        max_range = e
                new_SN_Array.append(max_range)
            elif selection == 'max_coverage_choose_uv':
                max_range = events[0]
                for e in events:
                    if (e.source == 'swift_uv' or e.source == 'uv' 
                        and e.filename not in bad_files):
                        max_range = e
                    elif (e.minwave) < (max_range.minwave):
                        max_range = e
                new_SN_Array.append(max_range)
            elif selection == 'choose_bluest':
                bluest = events[0]
                for e in events:
                    if e.minwave < bluest.minwave:
                        bluest = e
                new_SN_Array.append(bluest)
            elif selection == 'max_snr':
                max_snr = events[0]
                for e in events:
                    if (e.SNR != None and max_snr.SNR != None 
                        and abs(e.SNR) > abs(max_snr.SNR)):
                        max_snr = e
                new_SN_Array.append(max_snr)
            elif selection == 'accurate_phase':
                #this should be smarter
                ac_phase = events[0]
                for e in events:
                    if abs(e.phase - p_avg) < abs(ac_phase.phase - p_avg):
                        ac_phase = e
                new_SN_Array.append(ac_phase)
            elif selection == 'max_coverage_splice':
                max_range = events[0]
                for e in events:
                    if (e.maxwave - e.minwave) > (max_range.maxwave - max_range.minwave):
                        max_range = e
                splice_specs = []
                cur_min = max_range.minwave
                cur_max = max_range.maxwave
                print max_range.filename, max_range.minwave, max_range.maxwave
                max_spec_range = (max_range.wavelength[np.where((max_range.wavelength > max_range.minwave) 
                                                              & (max_range.wavelength < max_range.maxwave))])
                cur_minwave = max_range.minwave
                cur_maxwave = max_range.maxwave
                for e in events:
                    if e != max_range:
                        spec_range = (e.wavelength[np.where((e.wavelength > e.minwave) 
                                                          & (e.wavelength < e.maxwave))])
                        olap = 2.*len(np.intersect1d(spec_range, max_spec_range))
                        if olap < 500.:
                            if e.minwave < cur_minwave:
                                cur_minwave = e.minwave
                            if e.maxwave < cur_maxwave:
                                cur_maxwave = e.maxwave
                            max_spec_range = (max_range.wavelength[np.where((max_range.wavelength >cur_minwave)
                                                                          & (max_range.wavelength < cur_maxwave))])
                            splice_specs.append(e)

                new_SN_Array.append(max_range)
                for spec in splice_specs:
                    new_SN_Array.append(spec)
            elif selection == 'max_coverage_add_uv':
                max_range = events[0]
                for e in events:
                    if (e.maxwave - e.minwave) > (max_range.maxwave - max_range.minwave):
                        max_range = e
                add_specs = []
                best_add=None
                for e in events:
                    if e != max_range:
                        if e.minwave < 2500. and e.minwave < max_range.minwave:
                            add_specs.append(e)
                if len(add_specs)>0:
                    best_add = add_specs[0]
                for e in add_specs:
                    if best_add != None and e.minwave < best_add.minwave:
                        best_add = e

                new_SN_Array.append(max_range)
                if best_add != None:
                    new_SN_Array.append(best_add)

        SN_Array = new_SN_Array

    print len(SN_Array), "valid SNe found"

    for SN in SN_Array:
        #assign more attributes
        SN.phase_array = np.array(SN.flux)
        SN.dm15_array  = np.array(SN.flux)
        SN.red_array   = np.array(SN.flux)
        SN.vel         = np.array(SN.flux)

        nan_bool_flux = np.isnan(SN.flux)
        non_nan_data = np.where(nan_bool_flux == False)
        nan_data = np.where(nan_bool_flux == True)
        SN.flux[nan_data]         = np.nan
        SN.ivar[nan_data]         = 0.
        SN.phase_array[nan_data]  = np.nan
        SN.dm15_array[nan_data]   = np.nan
        SN.red_array[nan_data]    = np.nan
        SN.vel[nan_data]          = np.nan
        if SN.phase != None:
            SN.phase_array[non_nan_data] = SN.phase
        else:
            SN.phase_array[non_nan_data] = np.nan
        if get_phot:
            if SN.dm15_source != None:
                SN.dm15_array[non_nan_data] = SN.dm15_source
            elif SN.dm15_from_fits != None:
                SN.dm15_array[non_nan_data] = SN.dm15_from_fits
            else:
                SN.dm15_array[non_nan_data] = np.nan
        if SN.redshift != None:
            SN.red_array[non_nan_data] = SN.redshift
        else:
            SN.red_array[non_nan_data] = np.nan
        if SN.velocity != None:
            SN.vel[non_nan_data] = SN.velocity
        else:
            SN.vel[non_nan_data] = np.nan

        non_nan_data = np.array(non_nan_data[0])
        if len(non_nan_data) > 0:
            SN.x1 = non_nan_data[0]
            SN.x2 = non_nan_data[-1]
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
    """Counts the number of spectra contributing to the composite at any given wavelength.
    
        Args:
            SN_Array: An array of supernova objects.

        Returns:
            An array containing the the number of spectra in SN_Array that have
            flux at a given wavelength element.
    """
    spec_per_bin = []
    
    for i in range(len(SN_Array[0].flux)):
        count = 0
        for SN in SN_Array:
            if SN.flux[i] != 0 and SN.ivar[i] != 0:   
                count += 1
        spec_per_bin.append(count)
    
    return spec_per_bin
            
            
def optimize_scales(SN_Array, template, initial, scale_range=False, wave1=3000, wave2=6000):
    """Scales each unique supernova in SN_Array by minimizing the square residuals 
    between the supernova flux and the template flux. This also works for bootstrap 
    arrays (can contain repeated spectra) because the objects in SN_Array are not copies.

        Args:
            SN_Array: An array of supernova objects.
            template: The SN object that others will be scaled to
            initial: If True minimmize square residuals, if False minimize 
                square residuals weighted by inverse variance

        Returns:
            SN_Array: An array of supernova objects with flux, ivar, low_conf, 
                and up_conf scaled to the template
            scales: A list of the scales used corresponding to each SN object
    """
    scales = []
    unique_arr = list(set(SN_Array))
    guess = 1.0
    for uSN in unique_arr:
        guess = np.average(template.flux[template.x1:template.x2])/np.average(uSN.flux[uSN.x1:uSN.x2])
        u = opt.minimize(sq_residuals, guess, args = (uSN, template, initial, scale_range, wave1, wave2), 
                         method = 'Nelder-Mead').x
        scales.append(u)
        
    for i in range(len(unique_arr)):
        unique_arr[i].flux = scales[i]*unique_arr[i].flux
        unique_arr[i].ivar /= (scales[i])**2
        if len(unique_arr[i].low_conf) > 0 and len(unique_arr[i].up_conf) > 0:
            unique_arr[i].low_conf = scales[i]*unique_arr[i].low_conf
            unique_arr[i].up_conf = scales[i]*unique_arr[i].up_conf

    return SN_Array, scales
    
    
def sq_residuals(s,SN,comp, initial, scale_range, wave1, wave2):
    """Calculates the sum of the square residuals or weighted square residuals 
    between two supernova flux arrays using a given scale s.

    Args:
        s: The scaling factor
        SN: A supernova object 
        comp: Another supernova object (typically the composite spectrum)
        initial: If True minimmize square residuals, if False minimize 
            square residuals weighted by inverse variance.

    Returns:
        The sum of the square resiudal or weighted square residuals in the 
        region where the two spectra overlap.
    """
    s = np.absolute(s)
    #This is dumb but its works
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
        pos1 = 0
        pos2 = 0
    temp_flux = s*SN.flux

    if scale_range:
        wrange = np.where((SN.wavelength >= wave1) & (SN.wavelength < wave2))[0]
        pos1 = wrange[0]
        pos2 = wrange[-1]

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
    (creating the composite or bootstrapping). Applies masks the these data for 
    consistency and zero compensation. Returns the masks and data structures.
    
    Args:
        SN_Array: An array of supernova objects
        boot: If True, does not create data structures for redshift, phase,
            velocity, or dm15. 

    Returns:
        Masked 2D arrays containing spectral and photometric data for each SN object 
        to be combined as a function of wavelength.
    """

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
        fluxes.append(SN.flux)
        ivars.append(SN.ivar)
        if not boot:
            reds.append(SN.red_array)
            phases.append(SN.phase)
            ages.append(SN.phase_array)
            vels.append(SN.vel)
            dm15s.append(SN.dm15_array)
    
    fluxes = np.ma.masked_array(fluxes,np.isnan(fluxes))
    reds   = np.ma.masked_array(reds,np.isnan(reds))
    phases = np.ma.masked_array(phases,np.isnan(phases))
    ages   = np.ma.masked_array(ages,np.isnan(ages))
    vels   = np.ma.masked_array(vels,np.isnan(vels))
    dm15s  = np.ma.masked_array(dm15s,np.isnan(dm15s))
    dm15_ivars = np.ma.masked_array(dm15_ivars,np.isnan(dm15_ivars))
    red_ivars  = np.ma.masked_array(red_ivars,np.isnan(red_ivars))
    dm15_mask  = np.ma.masked_array(dm15_mask,np.isnan(dm15_mask))
    red_mask   = np.ma.masked_array(red_mask,np.isnan(red_mask))
    flux_mask = []
    ivar_mask = []
    
    return (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, 
            dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
                

def average(SN_Array, template, medmean, boot, fluxes, ivars, dm15_ivars, 
            red_ivars, reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, 
            dm15_mask, red_mask, find_RMSE=False):
    """Modifies the template supernova to be the inverse variance weighted 
    average of the scaled data.

        Args:
            SN_Array: An array of supernova objects
            template: The supernova object containing the composite spectrum
            medmean: If 1, do an inverse variance weighted average as a function 
                of wavelength. If 2, do a median.
            boot: If true, only combine data from flux arrays.
            *other Args*: The masked arrays generated from the 'mask' function.

        Keyword Arguments:
            find_RMSE: Determine the RMSE associated with the composite spectrum 
                and assign it as an attribute.

        Returns:
            A supernova object that has had it wavelength dependent attributes 
            combined.

    """
    temp_fluxes = []
    if medmean == 1: 
        template.flux  = np.ma.average(fluxes, weights=ivars, axis=0).filled(np.nan)
        if not boot:
            template.phase_array   = np.ma.average(ages, weights=ivars, axis=0).filled(np.nan)
            template.vel   = np.ma.average(vels, weights=ivars, axis=0).filled(np.nan)
            template.dm15_array  = np.ma.average(dm15s, weights=ivars, axis=0).filled(np.nan)
            template.red_array = np.ma.average(reds, weights = ivars, axis=0).filled(np.nan)
    if medmean == 2:
        template.flux = np.ma.median(fluxes, axis=0).filled(np.nan)
        if not boot:
            template.phase_array   = np.ma.median(ages, axis=0).filled(np.nan)
            template.vel   = np.ma.median(vels, axis=0).filled(np.nan)
            template.dm15_array  = np.ma.median(dm15s, axis=0).filled(np.nan)
            template.red_array = np.ma.median(reds, axis=0).filled(np.nan)
    
    #finds and stores the variance data of the template
    no_data   = np.where(np.sum(ivars, axis = 0)==0)
    template.ivar = np.sum(ivars, axis=0)
    template.ivar[no_data] = 0
    template.name = "Composite Spectrum"

    template.RMSE = []
    if find_RMSE:
        sq_diffs = np.ma.sum(ivars*(fluxes - template.flux)**2., axis=0).filled(np.nan)
        rmse = np.sqrt(sq_diffs/template.ivar)
        template.RMSE = rmse

    return template

def bootstrapping (SN_Array, samples, og_template, iters, medmean):
    """Creates a matrix of random sets of supernovae from the original sample 
    with the same size as the original sample. The number of samples is defined 
    by the user. Then creates the composite spectrum for each of these sets. 
    These data are used to contruct a confidence interval for the original sample. 

        Args:
            SN_Array: An array of supernova objects
            samples: The number of desired bootstrapping samples
            og_template: The supernova object corresponding to the original 
                composite spectrum.
            iters: The number of scaling iterations to make (this should be set
                to 1 if scaling technique is minimizing square residuals).
            medmean: If 1, do an inverse variance weighted average as a function 
                of wavelength. If 2, do a median.

        Returns:
            An array corresponding to the lower confidence interval, upper
            confidence interval, and an array of composite spectrum supernova 
            objects for every resampling of the data.

    """
    strap_matrix = np.random.random_sample((samples, len(SN_Array)))
    strap_matrix *= len(SN_Array)
    strap_matrix = strap_matrix.astype(int)   
    boot_arr = []
    boots = []
    boot = True
    
    cpy_array = []
    for SN in SN_Array:
        cpy_array.append(copy.deepcopy(SN))
        
    
    for i in range(len(strap_matrix)):
        boot_arr.append([])
        for j in range(len(strap_matrix[i])):
            boot_arr[i].append(copy.deepcopy(cpy_array[strap_matrix[i,j]]))

    for p in range(len(boot_arr)):
        # print p
        lengths = []
        for SN in boot_arr[p]:  
            lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
        boot_temp = [SN for SN in boot_arr[p] 
                     if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]
        boot_temp = copy.deepcopy(boot_temp[0])
        

        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
         flux_mask, ivar_mask, dm15_mask, red_mask) = mask(boot_arr[p], boot)
        for x in range(iters):
            new_SN_Array, scales = optimize_scales(boot_arr[p], boot_temp, True)
            (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
             flux_mask, ivar_mask, dm15_mask, red_mask) = mask(boot_arr[p], boot)
            template = average(boot_arr[p], boot_temp, medmean, boot, fluxes, 
                                ivars, dm15_ivars, red_ivars, reds, phases, ages, 
                                vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)
        boots.append(copy.deepcopy(template))

    print "scaling boots..."
    temp1, scales = optimize_scales(boots, og_template, True)

    #examine bootstrap samples
    # print "plotting..."
    # for SN in boots:
    #     plt.plot(SN.wavelength, SN.flux, 'g')
    # plt.plot(og_template.wavelength,og_template.flux, 'k', linewidth = 4)
    plt.show()
    
    print "computing confidence intervals..."
    resid = []
    percentile = erf(1/np.sqrt(2.))
    low_pc = 0.5 - percentile*0.5
    up_pc = 0.5 + percentile*0.5

    for SN in boots:
        resid.append(SN.flux - og_template.flux)
        
    resid_trans = np.transpose(resid)
    resid_sort = np.sort(resid_trans)
    arr = []
    new_resid = []
    for i in range(len(resid_sort)):
        if True in np.isfinite(resid_sort[i]):
            new_resid.append(resid_sort[i][np.isfinite(resid_sort[i])])
        else:
            new_resid.append(resid_sort[i])

    #plot a histogram of the residuals
    # for elem in resid_sort:
    #   for i in range (len(elem)):
    #       if np.isfinite(elem[i]):
    #           arr.append(elem[i])
    # plt.hist(arr,100)
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
    
    return np.asarray(low_arr), np.asarray(up_arr), boots        
    
def is_bad_data(SN, bad_files, bad_ivars):
    """Determines if a given supernova object has been flagged as having poor 
    data quality.

    Args:
        SN: A supernova object with an associated filename and event name
        bad_files: The list of bad files in questionable_spectra.py found to have
            poor data.
        bad_ivars: The list of files found to have incorrectly generated ivar spectra

    Returns:
        True if the data has been flagged. False if the data are good.

    TODO: Revisit the flagged data.
    """
    #2008ia is not deredshifted
    #2006X, 2002bf have anomalous velocity
    #1991bg variance spectrum seem wrong
    bad_sns = ['2002bf', '2006x', '1991bg', '2008ia']
    for el in bad_files:
        if SN.filename == el:
            return True
    for el in bad_ivars:
        if SN.filename == el:
            return True
    for el in bad_sns:
        if SN.name == el:
            return True
    return False

def remove_peculiars(SN_Array, file):
    """Removes objects from SN_Array that have been flagged as peculiar.

        Args:
            SN_Array: An array of supernova objects.
            file: A text file with a list of peculiar events.

        Returns:
            A new array of supernova objects without peculiar events
    """
    SN_Array_no_pecs = []
    count = 1
    with open(file) as f:
        names = np.loadtxt(f, dtype = str)
        for SN in SN_Array:
            if SN.name not in names:
                SN_Array_no_pecs.append(SN)
            else:
                # print count, SN.name
                count += 1

    return SN_Array_no_pecs

def check_host_corrections(SN_Array):
    """Removes objects from SN_Array that do not host extinction estimates.

        Args:
            SN_Array: An array of supernova objects.

        Returns:
            A new array of supernova objects that have host extinction estimates.
    """
    has_host_corr = []
    for SN in SN_Array:
        if SN.av_25 != None:
            has_host_corr.append(SN)
        elif SN.av_mlcs31 != None:
            has_host_corr.append(SN)
        elif SN.av_mlcs17 != None:
            has_host_corr.append(SN)
        # else:
        #   print SN.name
    SN_Array = has_host_corr
    print len(SN_Array), 'spectra with host corrections'
    return SN_Array

def apply_host_corrections(SN_Array, lengths, r_v = 2.5, verbose=True, low_av_test=None):
    """Correct supernova spectra in SN_Array for host extinction.

        Args:
            SN_Array: An array of supernova objects.
            lengths: An empty array to contain the length of each spectrum.

        Keyword Arguments:
            r_v: The value of R_V used for the light curve fit
            verbose: If True, print metadata of each supernova object as 
                they are corrected.

        Returns:
            An array of supernova objects whose flus and ivar have been corrected
            for host extinction.
    """
    corrected_SNs = []
    for SN in SN_Array:
        if verbose:
            print SN.name, SN.filename, SN.SNR, SN.dm15_source, SN.dm15_from_fits, \
                    SN.phase, SN.redshift, SN.source, SN.wavelength[SN.x1], \
                    SN.wavelength[SN.x2], SN.ned_host, SN.av_25, SN.minwave, SN.maxwave

        if SN.av_25 != None:
            if SN.av_25 > low_av_test or low_av_test == None:
                pre_scale = (1.e-15/np.average(SN.flux[SN.x1:SN.x2]))
                SN.flux = pre_scale*SN.flux
                SN.ivar = SN.ivar/(pre_scale*pre_scale)
                old_wave = SN.wavelength*u.Angstrom
                old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
                spec1d = Spectrum1D.from_array(old_wave, old_flux)
                old_ivar = SN.ivar*u.Unit('W m-2 angstrom-1 sr-1')
                new_flux, new_ivar = test_dered.host_correction(SN.av_25, r_v, SN.name, 
                                                                old_wave, old_flux, old_ivar)
                SN.flux = new_flux.value
                SN.ivar = new_ivar.value
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
            else:
                print SN.filename, 'has low reddening!'
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

        elif SN.av_mlcs31 != None:
            if SN.av_mlcs31 > low_av_test or low_av_test == None:
                pre_scale = (1.e-15/np.average(SN.flux[SN.x1:SN.x2]))
                SN.flux = pre_scale*SN.flux
                SN.ivar = SN.ivar/(pre_scale*pre_scale)
                old_wave = SN.wavelength*u.Angstrom
                old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
                spec1d = Spectrum1D.from_array(old_wave, old_flux)
                old_ivar = SN.ivar*u.Unit('W m-2 angstrom-1 sr-1')
                new_flux, new_ivar = test_dered.host_correction(SN.av_mlcs31, 3.1, SN.name, 
                                                                old_wave, old_flux, old_ivar)
                SN.flux = new_flux.value
                SN.ivar = new_ivar.value
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
            else:
                print SN.filename, 'has low reddening!'
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

        elif SN.av_mlcs17 != None:
            if SN.av_mlcs17 > low_av_test or low_av_test == None:
                pre_scale = (1.e-15/np.average(SN.flux[SN.x1:SN.x2]))
                SN.flux = pre_scale*SN.flux
                SN.ivar = SN.ivar/(pre_scale*pre_scale)
                old_wave = SN.wavelength*u.Angstrom
                old_flux = SN.flux*u.Unit('W m-2 angstrom-1 sr-1')
                spec1d = Spectrum1D.from_array(old_wave, old_flux)
                old_ivar = SN.ivar*u.Unit('W m-2 angstrom-1 sr-1')
                new_flux, new_ivar = test_dered.host_correction(SN.av_mlcs17, 1.7, SN.name, 
                                                                old_wave, old_flux, old_ivar)
                SN.flux = new_flux.value
                SN.ivar = new_ivar.value
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))
            else:
                print SN.filename, 'has low reddening!'
                corrected_SNs.append(SN)
                lengths.append(len(SN.flux[np.where(SN.flux != 0)]))

    SN_Array = corrected_SNs
    print len(SN_Array), 'SNs with host corrections'
    return SN_Array


def remove_tell_files(SN_Array):
    """Removes supernova objects flagged as having telluric contamination from 
    an array of supernova objects.

        Args:
            SN_Array: An array of supernova objects

        Returns:
            An array of supernova objects that have not been flagged for telluric
            contamination.
    """
    tell_files = tspec.tel_spec()
    has_tell_file = False 
    for SN in SN_Array:
        if SN.filename in tell_files:
            has_tell_file = True 

    if has_tell_file:
        SN_Array_wo_tell = []
        for SN in SN_Array:
            if SN.filename not in tell_files:
                SN_Array_wo_tell.append(copy.deepcopy(SN))
        return SN_Array_wo_tell
    else:
        return SN_Array

    

def create_composite(SN_Array, boot, template, medmean, gini_balance=False, aggro=.5):
    """Given an array of supernova objects, creates a composite spectrum and 
    estimate the 1-sigma bootstrap resampling error (if desired).

        Args:
            SN_Array: An array of supernova objects
            boot: If True, estimates the error via bootstrap resampling
            template: A supernova object to contain the composite spectrum 
                usually initialized as the spectrum with the largest wavelength
                range.
            medmean: If 1, do an inverse variance weighted average as a function 
                of wavelength. If 2, do a median.

        Keyword Arguments:
            gini_balance: If True, finds region of the composite spectrum dominated 
                by high SNR spectra, then deweights these spectra to have the the 
                median weight in that wavelength range. TODO: improve this 
                algorithm (currently produced representative composites, but
                deweighting straight to median might not be ideal). If False, ivar
                spectra are left as is. 

        Returns:
            A supernova object with the composite spectrum, and an array of 
            supernova objects that with composite spectra generated from the 
            bootstrap resampling step.
    """
    i = 0
    scales  = []
    iters = 1
    iters_comp = 1
    # boot = False

    if boot is True:
        bootstrap = 'y'
    else:
        bootstrap = 'n'
    print "Creating composite..."

    optimize_scales(SN_Array, template, True)

    # for SN in SN_Array:
    #   SN.ivar = np.ones(len(SN_Array[0].ivar))

    (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
     flux_mask, ivar_mask, dm15_mask, red_mask) = mask(SN_Array, False)

    #deweight high SNR spectra using gini coefficients
    if gini_balance:
        imbalanced = True
        gini_coeffs, num_specs, gini_ranges = gini.gini_coeffs(SN_Array)
        gini_range_meds = []
        print 'Gini balancing...'
        i=0
        first_iter = True
        prev_swaps = []
        while imbalanced:
            # print gini_coeffs
            # print num_specs
            gini_range_meds = []
            for g in gini_ranges:
                sn_meds = []
                for SN in SN_Array:
                    g_locs = np.where((SN.wavelength >= g[0]) & 
                                        (SN.wavelength < g[1]) & (SN.ivar > 0.))[0]
                    if np.nansum(SN.ivar[g_locs]) > 0.:
                        sn_meds.append(np.nansum(SN.ivar[g_locs]))
                # gini_range_meds.append(np.nanmedian(sn_meds)) #NOT IDEAL BUT ALWAYS CONVERGES
                if len(sn_meds) > 0:
                    gini_range_meds.append(aggro*np.amax(sn_meds)) #THIS DETERMINES HOW TO DEWEIGHT AND SHOULD BE IMPROVED
                else:
                    gini_range_meds.append(np.nan)
            # print gini_range_meds
            # (deweight_SNs, scale_dict, scale_ref_dict = 
            #     gini.calc_deweight_ranges(SN_Array, gini_coeffs, gini_ranges, 
            #                                 gini_range_meds, tol=.85))
            deweight_SNs, scale_dict, scale_ref_dict = \
            gini.calc_deweight_ranges(SN_Array, gini_coeffs, gini_ranges, 
                                        gini_range_meds, tol=.6)
            if len(deweight_SNs) == 0:
                imbalanced = False
            else:
                gini.deweight_biasing_SNe(deweight_SNs, scale_dict, scale_ref_dict)
                gini_coeffs, num_specs, gini_ranges = gini.gini_coeffs(SN_Array)
            i+=1
            # if i == 10:
            #   imbalanced = False
            first_iter = False
        print 'Balanced after', i, 'iterations'

    # qdb.plot_comp_and_all_spectra(template, SN_Array, show_ivar=True)
    for i in range(iters_comp):
        # SN_Array, scales = optimize_scales(SN_Array, template, False)
        SN_Array, scales = optimize_scales(SN_Array, template, True)
        # for SN in SN_Array:
        #   plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
        # plt.plot(template.wavelength, template.flux, 'k-', linewidth=4)
        # plt.show()
        # for SN in SN_Array:
        #   SN.ivar = np.ones(len(SN_Array[0].ivar))

        (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
         flux_mask, ivar_mask, dm15_mask, red_mask) = mask(SN_Array, False)
        template = average(SN_Array, template, medmean, False, fluxes, ivars, 
                            dm15_ivars, red_ivars, reds, phases, ages, vels, 
                            dm15s, flux_mask, ivar_mask, dm15_mask, red_mask, find_RMSE=True)
        # for SN in SN_Array:
        #   plt.plot(SN.wavelength, SN.flux)
        # plt.plot(template.wavelength, template.flux, 'k-', linewidth=4)
        # plt.show()
    print "Done."
    boots = None
    norm = 1./np.amax(template.flux[template.x1:template.x2])
    template.flux = template.flux*norm
    template.ivar = template.ivar/(norm**2.)
    if template.RMSE != None:
        template.RMSE = template.RMSE*(norm)

    for SN in SN_Array:
        SN.flux = SN.flux*norm
        SN.ivar = SN.ivar/(norm**2.)

        if SN.name == '1991bg':
            print 'Original:', np.nanmedian(SN.ivar[SN.x1:SN.x2])

    #create bootstrap composites
    if bootstrap is 'y':
        scales  = []
        print "Bootstrapping"
        samples = 100
        template.low_conf, template.up_conf, boots = \
            bootstrapping(SN_Array, samples, template, iters, medmean)
        up_diff = template.up_conf - template.flux
        low_diff = template.flux - template.low_conf

    nan_bool_flux = np.isnan(template.flux)
    non_nan_data = np.where(nan_bool_flux == False)
    non_nan_data = np.array(non_nan_data[0])
    if len(non_nan_data) > 0:
        template.x1 = non_nan_data[0]
        template.x2 = non_nan_data[-1]
        template.x2 += 1

    return template, boots
    
def main(Full_query, boot = 'nb', medmean = 1, make_corr=True, multi_epoch=False, 
        selection = 'max_coverage', gini_balance=False, aggro=.5, verbose=True, 
        low_av_test = None):
    """Main function. Finds spectra that satisfy the users query and creates a 
    composite spectrum based on the given arguments.
        
        Args:
            Full_query: The SQL query string 

        Keyword Arguments:
            boot: If 'b', estimate error via bootstrap resampling. If 'nb', create
                composite spectrum without bootstrapping (faster).
            medmean: If 1, do an inverse variance weighted average as a function 
                of wavelength. If 2, do a median.
            make_corr: If True, remove spectra that have been marked as 
                'questionable', peculiar events (Iax), and spectra that do not have 
                host extinction estimates.
            multi_epoch: If True, include all spectra for a given SN that satisify 
                the query. If False, choose one 1 spectrum per SN based on the 
                selection keyword.
            selection: If multi_epoch is False, this string defines the selection
                criteria for choosing a single spectrum from a SN. Options are:
                'maximum_coverage'(default): largest wavelength range
                'maximum_coverage_choose_uv': largest wavelength range but prioritize
                    hst and swift spectra
                'choose_bluest': smallest minimum wavelength
                'max_snr': highest signal to noise
                'accurate_phase': closest to middle of the phase bin. TODO: implement
                    this without parsing the query (currently requires uncommenting
                    code)
                'max_coverage_splice': allows multiple spectra from the same SN as 
                    long as overlap is < 500 A
            gini_balance: If True, finds region of the composite spectrum dominated 
                by high SNR spectra, then deweights these spectra to have the the 
                median weight in that wavelength range. TODO: improve this 
                algorithm (currently produced representative composites, but
                deweighting straight to median might not be ideal). If False, ivar
                spectra are left as is. 
            verbose: If True, print metadata of each supernova object in SN_Array.

        Returns:
            The supernova object containing the composite spectrum, the (potentially
            modified) SN_Array used to create the composite spectrum, and a list 
            supernova objects containing composite spectra generated from the 
            bootstrap resampling process.
    """
    SN_Array = []

    #Accept SQL query as input and then grab from the database
    print "SQL Query:", Full_query
    SN_Array = grab(Full_query, make_corr=make_corr, multi_epoch=multi_epoch, 
                    selection = selection)

    lengths = []
    SN_Array_wo_tell = remove_tell_files(SN_Array)
    print len(SN_Array) - len(SN_Array_wo_tell), 'spectra may have telluric contamination'

    SN_Array = apply_host_corrections(SN_Array, lengths, verbose=verbose, low_av_test=low_av_test)

    temp = [SN for SN in SN_Array if len(SN.flux[np.where(SN.flux!=0)]) == max(lengths)]

    try:
        composite = temp[0]
    except IndexError:
        print "No spectra found"
        exit()

    # spec_bin = spectra_per_bin(SN_Array)

    #finds range of useable data
    template = supernova()
    template = copy.deepcopy(composite)
    template.spec_bin = spectra_per_bin(SN_Array)

    #creates main composite
    i = 0
    scales  = []
    iters = 3
    iters_comp = 3
    if boot == 'b':
        boot = True
    else:
        boot = False
    # plt.figure(num = 2, dpi = 100, figsize = [30, 20], facecolor = 'w')
    # for i in range(len(SN_Array)):
    #     plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
    # plt.plot(template.wavelength, template.flux, 'k', linewidth = 4)
    # plt.show()

    #for updating one spectrum at a time
    # num_plots = len(SN_Array)
    # for i in range(num_plots):
    #   sub_sns = copy.copy(SN_Array[0:i+1])
    #   print SN_Array[i].filename, SN_Array[i].phase, SN_Array[i].source, SN_Array[i].SNR, np.average(SN_Array[i].ivar[SN_Array[i].x1:SN_Array[i].x2])
    #   optimize_scales(sub_sns, template, True)
    #   (fluxes, ivars, dm15_ivars, red_ivars, reds, phases, ages, vels, dm15s, 
    #    flux_mask, ivar_mask, dm15_mask, red_mask) = mask(sub_sns, boot)
    #   for j in range(iters_comp):
    #       SN_Array, scales = optimize_scales(SN_Array, template, False)
    #       template = average(sub_sns, template, medmean, boot, fluxes, ivars, dm15_ivars, red_ivars, 
    #                          reds, phases, ages, vels, dm15s, flux_mask, ivar_mask, dm15_mask, red_mask)

    #   plt.figure(num = 2, dpi = 100, figsize = [30, 15], facecolor = 'w')
    #   plt.subplot(2,1,1)
    #   # plt.plot(sub_sns[-1].wavelength, sub_sns[-1].flux)
    #   for SN in sub_sns:
    #       plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
    #   plt.plot(template.wavelength[SN.x1:SN.x2], template.flux[SN.x1:SN.x2], 'k', linewidth = 4)
    #   plt.subplot(2,1,2)
    #   # plt.plot(sub_sns[-1].wavelength, sub_sns[-1].ivar)
    #   for SN in sub_sns:
    #       plt.plot(SN.wavelength[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2])
    #       r = sa.measure_si_ratio(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
    #       print 'SN Si Ratio: ', r
    #   r = sa.measure_si_ratio(template.wavelength[template.x1:template.x2], template.flux[template.x1:template.x2], vexp = .001)
    #   print 'Comp Si Ratio: ', r
    #   plt.show()


    template, boots = create_composite(SN_Array, boot, template, 
                                        medmean, gini_balance=gini_balance, aggro=aggro)

    return template, SN_Array, boots

if __name__ == "__main__":
    main()
