import os
import glob
from specutils import extinction as ex
#import astroquery
#from astroquery.ned import Ned
#from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import math
from datafidelity import *  # Get variance from the datafidelity outcome

mn.patch()

#Sets up some lists for later
SN_Array = []
full_array = []

def grab(sql_input, Full_query):
    print "Collecting data..."
    SN_Array = []
    cur.execute(sql_input)
    for row in cur:
        if 1 == 1:
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

sql_input = 'Stuff'
X = []
Y = []

SN_Array = grab(sql_input, sql_input)

for SN in SN_Array:
    X.append(SN.B_minus_v)
    Y.append(SN.velocity)

plt.plot(X, Y)
plt.show()
