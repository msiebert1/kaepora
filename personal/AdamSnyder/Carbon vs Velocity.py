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

    return indices

# Connect to database
con = sq3.connect('..\..\data\SNe.db')
cur = con.cursor()
sql_input = 'SELECT * FROM Supernovae WHERE Velocity > -20 AND Phase BETWEEN -5 AND 5 AND Dm15 BETWEEN 1.0 AND 1.5'
X = []
Y = []
carbon = []
no_carbon = []
maybe = []

# Get SN data from databae
SN_Array_old = grab(sql_input, sql_input)
SN_Array = np.array([])

# Take only unique SN for use in plots
unique_SN = findunique(SN_Array_old)
for i in range(len(unique_SN)):
    if SN_Array_old[unique_SN[i]].carbon == None:
        continue
    else:
        SN_Array = np.append(SN_Array_old[unique_SN[i]], SN_Array)

for SN in SN_Array:
    if SN.carbon == "A":
        carbon.append(SN.velocity)
    elif SN.carbon == "N":
        no_carbon.append(SN.velocity)
    else:
        maybe.append(SN.velocity)

print sum(SN.velocity <= -12 for SN in SN_Array)
print sum(SN.velocity >-12 for SN in SN_Array)

bins = np.linspace(-16, -9, 20)

# Create figure
plt.title('Carbon Abundance vs. Si II Velocity')
plt.ylabel('Frequency')
plt.xlabel(r'$v_\mathrm{abs}$ (Si II $\lambda 6355$) [$10^3$ km/s]')

plt.gca().invert_xaxis()
plt.hist(no_carbon, bins, color = "Blue", alpha =0.5, label = "No Carbon")
plt.hist(maybe, bins, color = "Green", alpha = 0.5, label = "Likely Carbon")
plt.hist(carbon, bins, color = "Red", alpha = 0.5, label = "Carbon")
plt.legend()
plt.savefig('carbon-velocity-scatter.png')

#plt.scatter(X, Y)
plt.show()
