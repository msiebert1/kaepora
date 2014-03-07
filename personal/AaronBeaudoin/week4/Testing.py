import sqlite3 as sq3
import matplotlib as plt
import numpy as np
import msgpack as msg
import msgpack_numpy as mn
from io import BytesIO

mn.patch()


class supernova(object):
    """Attributes can be added"""

con = sq3.connect('../../../../../SNe.db')
cur = con.cursor()
SN_Array = []
compare_spectrum = []
max_light = []
max_light = np.loadtxt("MaxSpectra.dat", dtype = 'str', delimiter = " ", skiprows = 1)
names = []
j = 1
#Reads the data of max light spectra so only one per supernova is used
for line in max_light[0:100]:
	SN = supernova()
	SN.address = line[1]
	SN.name = line[0]
	data = np.loadtxt(SN.address)
	SN.wavelength = data[:,0]
	SN.flux = data[:,1]
	SN.residual = np.zeros(len(data[:,1]))   
	SN.age = line[2]
	error = data[:,2]
	if not all(x == 0.0 for x in error):
		SN.error = error
	SN_Array.append(SN)
	#print j
	j += 1
print j, 'supernovae fluxed'
#Reads from database, may be changed later
print "Reading properties from database..."
j = 0
buf = BytesIO()
for SN in SN_Array[0:10]:
    for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave, Spectra FROM Supernovae'):
	if row[0] in SN.address:
	    SN.filename = row[0]
	    SN.redshifts = row[2]
	    SN.minwave = row[3]
	    SN.maxwave = row[4]

            print SN.spectrum
	    names.append(SN.name)
	    #print j
	    #print SN.redshifts
	    j += 1
	else:
	    continue