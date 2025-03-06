import numpy as np 
from astropy.io import fits
import matplotlib.pyplot as plt 
from math import *
from specutils.io import read_fits

sky = fits.open('../personal/AdamSnyder/kecksky.fits')
crval = sky[0].header['CRVAL1']
delta = sky[0].header['CDELT1']
pix = sky[0].header['CRPIX1']
print pix
skyflux = sky[0].data[0]
start = crval 
stop = crval + ceil(len(skyflux)*delta)
skywave = [(start+delta*i) for i in range(len(skyflux))]

plt.plot(skywave, skyflux)
plt.show()