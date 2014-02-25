import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage

SN=np.genfromtxt('sn2001fu-20011107.51-fast.flm')

wavelength = SN[:,0]
flux = SN[:,1]
variance = SN[:2]

nflux = flux

sflux = scipy.ndimage.filters.gaussian_filter(nflux,1)
rat = nflux/sflux

#plt.plot(rat)
 
for i in range(len(nflux)):
    if rat[i] > 1.05:
            nflux[i] = sflux[i]

plt.plot(wavelength, nflux)
plt.plot(wavelength, flux, 'k')
plt.show()