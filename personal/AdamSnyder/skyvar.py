import numpy as np
import pyfits
import math
import datafidelity as df
import matplotlib.pyplot as plt

SN = np.genfromtxt('sn1996C-19960217.48-fast.flm')

wavelength = np.array(SN[:,0])
flux = np.array(SN[:,1])
error = np.array(SN[:, 2])

ivar = df.genivar(wavelength, flux, error)

ivar2 = 1 / (sm_error**2)

plt.plot(wavelength, ivar)
plt.show()
plt.plot(wavelength, ivar2 )
plt.show()
