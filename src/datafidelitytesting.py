import datafidelity as df
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *

#SN = np.genfromtxt('../../data/spectra/bsnip/sn2004bz-20040613.453-ui.flm')
SN = np.genfromtxt('../personal/AdamSnyder/sn1996C-19960217.48-fast.flm')

wavelength = SN[:, 0]
flux = SN[:, 1]

sky = pyfits.open('../personal/AdamSnyder/kecksky.fits')
    
crval = sky[0].header['CRVAL1']
delta = sky[0].header['CDELT1']
skyflux = sky[0].data[0]
start = crval - ceil(0.5*len(skyflux)*delta)
stop = crval + ceil(0.5*len(skyflux)*delta)
skywave = [(start+delta*i) for i in range(len(skyflux))]

plt.plot(skywave, skyflux)

print wavelength[0], wavelength[-1]

plt.plot(wavelength, flux)

ivar = df.genivar(wavelength, flux)

ivar_new = df.clip(wavelength, flux, ivar)

plt.plot(wavelength, flux)
plt.show()
