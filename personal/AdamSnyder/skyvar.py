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
con = sq3.connect('../../data/SNe.db')
cur = con.cursor()

#SN = np.genfromtxt('../../data/spectra/bsnip/sn2008dx-20080707.217-ui.flm')

cur.execute('SELECT * FROM Supernovae WHERE Phase BETWEEN -3 AND 3 AND Velocity BETWEEN -12 AND -8')
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

SN_Array = [SN for SN in SN_Array if hasattr(SN, 'wavelength')]
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'ivar')]
SN_Array = [SN for SN in SN_Array if SN.phase != None]
SN_Array = [SN for SN in SN_Array if SN.redshift != None]

print SN_Array[0].filename

plt.plot(SN_Array[0].wavelength, SN_Array[0].flux)
plt.show()
