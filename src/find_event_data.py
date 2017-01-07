import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np

con = sq3.connect('../data/SNe_2.db')
cur = con.cursor()

# np.set_printoptions(threshold=np.nan)
mn.patch()

SN_Array = []
full_array = []
compare_spectrum = []

def find_data(event):
	event = "'" + event + "'"
	query = "SELECT * FROM Supernovae where SN = " + event
	print query
	SN_Array = grab_event_data(query)
	return SN_Array

def find_all_data():
    query = "SELECT * FROM Supernovae "
    SN_Array = grab_event_data(query)
    return SN_Array

class supernova(object):
    """Contains all spectral data provided by the associated file.
       Attributes can be added
    """

def grab_event_data(sql_input):
    """Pulls in all columns from the database for the selected query. 
       Replaces all NaN values with 0. Returns the array of supernova objects 
       with the newly added attributes.
    """
    print "Collecting data..."
    SN_Array = []
    multi_epoch = True

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
