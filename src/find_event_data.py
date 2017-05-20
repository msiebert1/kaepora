import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import magnitudes as mag

con = sq3.connect('../data/SNe_14_phot_1.db')
cur = con.cursor()

# np.set_printoptions(threshold=np.nan)
mn.patch()

SN_Array = []
full_array = []
compare_spectrum = []

def find_data(event):
    event = "'" + event + "'"
    query = "SELECT * FROM Supernovae where SN = " + event
    phot_query = "SELECT * FROM Photometry where SN = " + event
    print query
    print phot_query
    SN_Array = grab_event_data(query)
    PH_Array = grab_event_phot_data(phot_query)
    return SN_Array, PH_Array

def find_all_data():
    query = "SELECT * FROM Supernovae "
    SN_Array = grab_event_data(query)
    return SN_Array

def find_csp_other_data():
    query = "SELECT * FROM Supernovae where source = 'csp' or source = 'other'"
    SN_Array = grab_event_data(query)
    return SN_Array

class supernova(object):
    """Contains all spectral data provided by the associated file.
       Attributes can be added
    """

class photometry(object):
    """Contains all photometric data for the specified event
    """

def grab_event_phot_data(sql_input):

    PH_Array = []
    cur.execute(sql_input)
    for row in cur:
        PH           = photometry()
        PH.name      = row[0]
        PH.re        = row[1]
        PH.dec       = row[2]
        PH.zCMB_salt, PH.e_zCMB_salt, PH.Bmag_salt, PH.e_Bmag_salt, PH.s_salt, PH.e_s_salt, PH.c_salt, PH.e_c_salt, PH.mu_salt, PH.e_mu_salt = row[3:13]
        PH.zCMB_salt2, PH.e_zCMB_salt2, PH.Bmag_salt2, PH.e_Bmag_salt2, PH.x1_salt2, PH.e_x1_salt2, PH.c_salt2, PH.e_c_salt2, PH.mu_salt2, PH.e_mu_salt2 = row[13:23]
        PH.zCMB_mlcs31, PH.e_zCMB_mlcs31, PH.mu_mlcs31, PH.e_mu_mlcs31, PH.delta_mlcs31, PH.e_delta_mlcs31, PH.av_mlcs31, PH.e_av_mlcs31 = row[23:31]
        PH.zCMB_mlcs17, PH.e_zCMB_mlcs17, PH.mu_mlcs17, PH.e_mu_mlcs17, PH.delta_mlcs17, PH.e_delta_mlcs17, PH.av_mlcs17, PH.e_av_mlcs17 = row[31:39]
        PH.glon_host, PH.glat_host, PH.cz_host, PH.czLG_host, PH.czCMB_host, PH.mtype_host, PH.xpos_host, PH.ypos_host, PH.t1_host, PH.filt_host, PH.Ebv_host = row[39:50]
        PH.zCMB_lc, PH.zhel_lc, PH.mb_lc, PH.e_mb_lc, PH.c_lc, PH.e_c_lc, PH.x1_lc, PH.e_x1_lc, PH.logMst_lc, PH.e_logMst_lc, PH.tmax_lc, PH.e_tmax_lc, PH.cov_mb_s_lc, PH.cov_mb_c_lc, PH.cov_s_c_lc, PH.bias_lc = row[50:66]
        PH.av_25 = row[66]
        PH.dm15_cfa = row[67]
        PH.light_curves = msg.unpackb(row[68])

        PH_Array.append(PH)    

    return PH_Array


def grab_event_data(sql_input):
    """Pulls in all columns from the database for the selected query. 
       Replaces all NaN values with 0. Returns the array of supernova objects 
       with the newly added attributes.
    """
    print "Collecting data..."
    SN_Array = []
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
        phot         = row[17]
        SN.phot      = phot
        SN.mjd          = row[18]
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
            SN.x2 += 1
        else:
            SN.x1 = 0
            SN.x2 = 0
                    
    print "Arrays cleaned"
    return SN_Array

# arr = find_data('2005cf')
# # mags = mag.ab_mags(arr)
# # print mags
# for SN in arr:
#     print SN.filename, SN.mjd, type(SN.mjd)