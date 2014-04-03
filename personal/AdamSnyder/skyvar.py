import numpy as np
import matplotlib.pyplot as plt
import pyfits
import math
from scipy import interpolate
import sqlite3 as sq3

def skyvar(wave, sky = 'kecksky.fits'):

    if False:
        sky = pyfits.open("kecksky.fits")
        crval = sky[0].header['CRVAL1']
        delta = sky[0].header['CDELT1']
        sky_flux = sky[0].data[0]
        print "Keck Sky used"
        
    else:
        sky = pyfits.open("licksky.fits")
        crval = sky[0].header['CRVAL1']
        delta = sky[0].header['CDELT1']        
        sky_flux = sky[0].data
        print "Lick sky used"

    start = crval - math.ceil(0.5*len(sky_flux)*delta)
    stop = crval + math.ceil(0.5*len(sky_flux)*delta)

    sky_wave = [(start+delta*i) for i in range(len(sky_flux))]

    plt.plot(sky_wave, sky_flux)
    plt.show()

    return new_sky

con = sq3.connect('../../../First/SNe.db')
cur = con.cursor()
maxdif = 0.0
start = 0.0
stop = 0.0
lowSN = []
highSN = []

# Iterate over entire database to find the Supernovae with greatest wavelength
for row in cur.execute('SELECT * FROM Supernovae'):
    dif = row[6] - row[5]
    if dif > maxdif:
        maxdif = row[6] - row[5]
        start = row[5]
        stop = row[6]
        SN = row[0]

minwav = start
maxwav = stop
maxwavstart = 0.0

for row in cur.execute('SELECT * FROM Supernovae'):
    if row[5] < start:
        lowSN.append(row[0])
    if row[6] > stop:
        highSN.append(row[0])
    if row[5] < minwav:
        print row[0]
        minwav = row[5]
    if row[6] > maxwav:
        print row[0]
        maxwav = row[6]
        maxwavstart = row[5]
        
# Get flux data of supernovae

#flux = []
#wavelength = []

# Generate continuum

#residual = flux - continuum


print SN, maxdif, start, stop
print len(lowSN), len(highSN), minwav, maxwav, maxwavstart
    

