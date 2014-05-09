import os
#import glob
#from specutils import extinction as ex
#import astroquery
#from astroquery.ned import Ned
#from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import math
#from datafidelity import *  # Get variance from the datafidelity outcome
import msgpack as msg
import msgpack_numpy as mn
import sqlite3 as sq3

mn.patch()

class supernova():
    """stuff"""

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
    return SN_Array

def findunique(SN_Array):
    Names = []
    indices = []
    for i in range(len(SN_Array)):
        if SN_Array[i].name in Names:
            pass
        else:
            Names.append(SN_Array[i].name)
            indices.append(i)

    print indices

    return indices

# Connect to database
con = sq3.connect('..\..\data\SNe.db')
cur = con.cursor()
sql_input = 'SELECT * FROM Supernovae WHERE Morphology BETWEEN 1 AND 11 AND Velocity > -20 AND Phase BETWEEN -3 AND 3 AND Dm15 BETWEEN 1.0 AND 1.5'
X = []
Y = []
high_morph = []
low_morph = []
high_vel = []
low_vel = []

SN_Array_old = grab(sql_input, sql_input)
SN_Array = np.array([])

unique_SN = findunique(SN_Array_old)

for i in range(len(unique_SN)):
    SN_Array = np.append(SN_Array_old[unique_SN[i]], SN_Array)

print len(SN_Array)

for SN in SN_Array:
    Y.append(SN.morph)
    X.append(SN.velocity)
    if SN.velocity <= -12:
        high_vel.append(SN.velocity)
        high_morph.append(SN.morph)
    elif SN.velocity > -12:
        low_vel.append(SN.velocity)
        low_morph.append(SN.morph)

print len(high_vel), len(low_vel)

# Create figure
fig = plt.figure()
plt.title('Galaxy Morphology vs. Si II line Velocity')
plt.ylabel(r'$v_\mathrm{abs}$ (Si II $\lambda 6355$) [$10^3$ km/s]')
plt.xlabel(r'Host Galaxy Morphology')

# Create best fit line
#fit = np.polyfit(X,Y,1)
#fit_fn = np.poly1d(fit)
#fit_label = 'B-V = {0:.4f} {1:.4f} x (v / 10^3 km/s)'.format(fit[1], fit[0])

#calculate Pearson correlation coefficient
#correlation = np.corrcoef(X,Y)[0, 1]
#print correlation

plt.gca().invert_yaxis()
plt.plot(high_morph, high_vel, 'ro', low_morph, low_vel, 'bo')
plt.savefig('morph-velocity-scatter.png')

#plt.scatter(X, Y)
plt.show()
